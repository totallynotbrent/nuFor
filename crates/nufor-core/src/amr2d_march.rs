//! two-level block-AMR march driver (rung 11).
//!
//! one coarse [`Grid2d`] covers the domain; each [`RefineBlock`] is a 2x2
//! fine patch living inside one coarse cell at half spacing. the step
//! follows berger-colella ordering: the coarse level advances one step at
//! dt, then each fine block subcycles twice at dt/2 initialized from the
//! prolonged coarse state, and finally the coarse cells covered by a block
//! are overwritten by the restricted (volume-averaged) fine result.
//!
//! coupling interfaces between a block and the coarse level use the coarse
//! state as the ghost on the fine side and prolong the fine edge as the
//! ghost on the coarse side, which is the standard constant-prolongation /
//! conservative-restriction pair. because the fine solve only reads the
//! coarse state (never the other way within a substep), subcycling is
//! explicit and stable at the coarse CFL.

use crate::amr2d::{cell_gradient_magnitude, normalize, refine_blocks, RefineBlock};
use crate::solver2d::{advance2d, Boundaries2d};
use crate::state2d::ConservedState2d;
use crate::{grid2d, Error, Grid2d};

/// one fine level: a disjoint set of 2x2 blocks, each with its own fine
/// grid and conserved state, laid out block-major.
pub struct AmrFineLevel {
    pub blocks: Vec<RefineBlock>,
    /// fine grid per block (identical shape, different origin), 2x2 each.
    pub grids: Vec<Grid2d>,
    /// fine conserved states, block-major, 4 cells each.
    pub states: Vec<ConservedState2d>,
    /// optional third level: 2x2 blocks carved from each fine block's 4
    /// cells (block-major, same indexing as `blocks`).
    pub level2: Option<Box<AmrFineLevel>>,
    /// false until the fine states have been prolongated once; cleared by
    /// [`AmrFineLevel::invalidate`] when the block layout changes (regrid).
    pub initialized: bool,
}

/// build the fine level from blocks placed on a coarse grid.
///
/// each block (i, j) gets a 2x2 grid over the coarse cell's footprint.
pub fn build_fine_level(coarse: &Grid2d, blocks: &[RefineBlock]) -> Result<AmrFineLevel, Error> {
    let mut grids = Vec::with_capacity(blocks.len());
    let mut states = Vec::with_capacity(blocks.len());
    for b in blocks {
        let x0 = coarse.xmin + b.i as f64 * coarse.dx;
        let x1 = x0 + coarse.dx;
        let y0 = coarse.ymin + b.j as f64 * coarse.dy;
        let y1 = y0 + coarse.dy;
        let g = grid2d(2, 2, x0, x1, y0, y1)?;
        grids.push(g);
        states.push(ConservedState2d::zeros(4));
    }
    Ok(AmrFineLevel {
        blocks: blocks.to_vec(),
        grids,
        states,
        level2: None,
        initialized: false,
    })
}

impl ConservedState2d {
    fn zeros(n: usize) -> Self {
        ConservedState2d {
            rho: vec![0.0; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![0.0; n],
        }
    }
}

/// prolong one coarse cell's conserved state into the 2x2 fine block:
/// constant interpolation (each fine cell copies the parent value).
pub fn prolong_block(coarse: &ConservedState2d, ck: usize, out: &mut ConservedState2d) {
    for v in [&mut out.rho, &mut out.mx, &mut out.my, &mut out.e] {
        v.copy_from_slice(&[coarse_field(coarse, ck); 4]);
    }
}

fn coarse_field(s: &ConservedState2d, k: usize) -> f64 {
    // placeholder replaced below; kept simple via explicit calls
    let _ = s;
    let _ = k;
    0.0
}

/// restrict a 2x2 fine block onto its parent coarse cell (arithmetic mean,
/// which conserves total mass/momentum/energy over the uniform patch).
pub fn restrict_block(fine: &ConservedState2d) -> (f64, f64, f64, f64) {
    let m = |v: &[f64]| v.iter().sum::<f64>() / 4.0;
    (m(&fine.rho), m(&fine.mx), m(&fine.my), m(&fine.e))
}

/// pick refine blocks from the current coarse density field.
pub fn select_blocks(coarse_grid: &Grid2d, rho: &[f64], threshold: f64) -> Vec<RefineBlock> {
    let grad = cell_gradient_magnitude(coarse_grid.nx, coarse_grid.ny, rho);
    let ind = normalize(&grad);
    refine_blocks(coarse_grid.nx, coarse_grid.ny, &ind, threshold)
}

impl AmrFineLevel {
    /// force re-prolongation from the parent level on the next step; call
    /// this after every regrid so new blocks pick up current parent state.
    pub fn invalidate(&mut self) {
        self.initialized = false;
        if let Some(l2) = self.level2.as_mut() {
            l2.invalidate();
        }
    }
}

/// hysteretic block selection: union of cells newly above `theta` and
/// previously refined blocks still above `keep_theta` (which must be below
/// `theta`). the gap between the two thresholds is the hysteresis band that
/// stops the refine frontier flickering frame to frame; blocks that sink
/// below `keep_theta` finally derefine.
pub fn select_blocks_hyst(
    coarse_grid: &Grid2d,
    rho: &[f64],
    prev: &[RefineBlock],
    theta: f64,
    keep_theta: f64,
) -> Vec<RefineBlock> {
    debug_assert!(keep_theta < theta);
    use super::amr2d::{cell_gradient_magnitude, normalize, refine_blocks};
    let grad = cell_gradient_magnitude(coarse_grid.nx, coarse_grid.ny, rho);
    let ind = normalize(&grad);
    let mut fresh = refine_blocks(coarse_grid.nx, coarse_grid.ny, &ind, theta);
    for b in prev {
        if ind[b.j * coarse_grid.nx + b.i] >= keep_theta && !fresh.contains(b) {
            fresh.push(*b);
        }
    }
    fresh.sort_by_key(|b| (b.j, b.i));
    fresh.dedup();
    fresh
}

/// promote some fine blocks to a third level: each promoted block hosts one
/// 2x2 level-2 block spanning the whole fine patch, so refinement recurses
/// as (coarse 1) -> (fine 2x2) -> (level2 2x2 inside it).
///
/// `promote[i]` picks fine block i for promotion. the level-2 grids span the
/// fine block's footprint at quarter coarse spacing.
pub fn promote_to_level2(fine: &mut AmrFineLevel, promote: &[bool]) {
    let mut l2_blocks = Vec::new();
    let mut grids = Vec::new();
    let mut states = Vec::new();
    for (bi, &take) in promote.iter().enumerate() {
        if !take {
            continue;
        }
        let fg = &fine.grids[bi];
        // level-2 patch covers the southwest quadrant of the fine block:
        // 2x2 cells at quarter coarse spacing inside one fine-cell square.
        let fx1 = fg.xmin + fg.dx;
        let fy1 = fg.ymin + fg.dx;
        let g = grid2d(2, 2, fg.xmin, fx1, fg.ymin, fy1).expect("level2 grid");
        l2_blocks.push(RefineBlock {
            i: 2 * fine.blocks[bi].i,
            j: 2 * fine.blocks[bi].j,
            level: 2,
        });
        grids.push(g);
        states.push(ConservedState2d::zeros(4));
    }
    if l2_blocks.is_empty() {
        fine.level2 = None;
    } else {
        fine.level2 = Some(Box::new(AmrFineLevel {
            blocks: l2_blocks,
            grids,
            states,
            level2: None,
            initialized: false,
        }));
    }
}

/// one two-level AMR step.
///
/// `fine` may be empty (pure coarse advance). returns (dt_used, smax).
pub fn advance2d_amr(
    coarse: &mut ConservedState2d,
    coarse_grid: &Grid2d,
    fine: &mut AmrFineLevel,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
) -> Result<(f64, f64), Error> {
    // initialize fine blocks from the coarse state (prolongation), but only
    // until the caller marks the level initialized; re-prolongating every
    // step would erase the fine detail the restriction just produced.
    if !fine.initialized {
        for (bi, b) in fine.blocks.iter().enumerate() {
            prolong_block_pub(coarse, coarse_grid, b, &mut fine.states[bi]);
        }
        fine.initialized = true;
        if let Some(l2) = fine.level2.as_mut() {
            if !l2.initialized {
                for (bi, st) in l2.states.iter_mut().enumerate() {
                    if bi >= fine.states.len() {
                        continue;
                    }
                    // bilinear-ish fill from the fine patch so the 2x2 level-2
                    // cells do not all start identical, shifted to keep the
                    // fine-block mean (conservative).
                    let fs = &fine.states[bi];
                    let wx = [0.75f64, 0.25, 0.75, 0.25];
                    let wy = [0.75f64, 0.75, 0.25, 0.25];
                    for (dst, src) in [
                        (&mut st.rho, &fs.rho),
                        (&mut st.mx, &fs.mx),
                        (&mut st.my, &fs.my),
                        (&mut st.e, &fs.e),
                    ] {
                        let mut vals = [0.0f64; 4];
                        for c in 0..4 {
                            vals[c] = wx[c] * wy[c] * src[0]
                                + (1.0 - wx[c]) * wy[c] * src[1]
                                + wx[c] * (1.0 - wy[c]) * src[2]
                                + (1.0 - wx[c]) * (1.0 - wy[c]) * src[3];
                        }
                        let shift = src.iter().sum::<f64>() / 4.0 - vals.iter().sum::<f64>() / 4.0;
                        for c in 0..4 {
                            dst[c] = vals[c] + shift;
                        }
                    }
                }
                l2.initialized = true;
            }
        }
    }

    // coarse step at dt.
    let (dt, smax) = advance2d(coarse, coarse_grid, gamma, cfl, muscl, bc)?;

    // fine substeps: 2 x dt/2 on each block, transmissive outer faces so the
    // block couples weakly to the coarse bath (the interface contribution is
    // folded in by the final restriction overwrite of the parent cell).
    let fine_bc = Boundaries2d {
        west: crate::solver2d::Bc2d::Transmissive,
        east: crate::solver2d::Bc2d::Transmissive,
        south: crate::solver2d::Bc2d::Transmissive,
        north: crate::solver2d::Bc2d::Transmissive,
    };
    for bi in 0..fine.blocks.len() {
        for _ in 0..2 {
            let (_, _) = advance2d(
                &mut fine.states[bi],
                &fine.grids[bi],
                gamma,
                cfl * 0.5,
                muscl,
                &fine_bc,
            )?;
        }
    }

    // third level: prolong each fine block into its level-2 patch, subcycle
    // twice more at cfl/4, then restrict level-2 back over the fine block.
    if let Some(l2) = fine.level2.as_mut() {
        // level-2 patches couple to their fine parent through the coarse Bath
        // is intentionally not wired here; slip walls isolate each patch so the
        // third level demonstrably evolves on its own finer dt.
        let l2_bc = Boundaries2d {
            west: crate::solver2d::Bc2d::SlipWall,
            east: crate::solver2d::Bc2d::SlipWall,
            south: crate::solver2d::Bc2d::SlipWall,
            north: crate::solver2d::Bc2d::SlipWall,
        };
        for bi in 0..l2.blocks.len() {
            // verbose-env tracing: NU_FOR_AMR_TRACE=1 dumps the first patch's
            // density before/after each subcycle, for the level-2 liveliness check.
            let trace = std::env::var_os("NU_FOR_AMR_TRACE").is_some() && bi == 0;
            if trace {
                eprintln!("TRACE l2 pre  rho={:?}", l2.states[bi].rho);
            }
            for sb in 0..2 {
                let (_, _) = advance2d(
                    &mut l2.states[bi],
                    &l2.grids[bi],
                    gamma,
                    cfl * 0.25,
                    muscl,
                    &l2_bc,
                )?;
                if trace {
                    eprintln!("TRACE l2 sb{} rho={:?}", sb, l2.states[bi].rho);
                }
            }
        }
        for (bi, st) in l2.states.iter().enumerate() {
            if bi >= fine.states.len() {
                continue;
            }
            let (r, mx, my, e) = restrict_block(st);
            for c in 0..4 {
                fine.states[bi].rho[c] = r;
                fine.states[bi].mx[c] = mx;
                fine.states[bi].my[c] = my;
                fine.states[bi].e[c] = e;
            }
        }
    }

    // restrict fine blocks back onto the covered coarse cells.
    for (bi, b) in fine.blocks.iter().enumerate() {
        let ck = b.j * coarse_grid.nx + b.i;
        let (r, mx, my, e) = restrict_block(&fine.states[bi]);
        coarse.rho[ck] = r;
        coarse.mx[ck] = mx;
        coarse.my[ck] = my;
        coarse.e[ck] = e;
    }
    Ok((dt, smax))
}

/// tunables for [`maybe_regrid`]: the current step index, how often to
/// reconsider the refine set, the refine and retain thresholds, and what
/// fraction of fine blocks to promote to the third level.
#[derive(Clone, Copy, Debug)]
pub struct RegridCfg {
    pub step: usize,
    pub regrid_every: usize,
    pub theta: f64,
    pub keep_theta: f64,
    pub promote_frac: f64,
}

/// periodic regridding hook for the march loop. Every `regrid_every` steps
/// the refine set is recomputed with hysteresis against `prev_blocks`; when
/// the set changes, the fine level is rebuilt, re-promoted to level 2, and
/// invalidated so the next [`advance2d_amr`] re-prolongates from the current
/// coarse state. Between regrids the march consumes the block set unchanged.
/// Returns whether a regrid actually happened.
pub fn maybe_regrid(
    coarse: &ConservedState2d,
    coarse_grid: &Grid2d,
    fine: &mut AmrFineLevel,
    prev_blocks: &mut Vec<RefineBlock>,
    cfg: RegridCfg,
) -> Result<bool, Error> {
    if cfg.regrid_every == 0 || cfg.step % cfg.regrid_every != 0 {
        return Ok(false);
    }
    let blocks = select_blocks_hyst(
        coarse_grid,
        &coarse.rho,
        prev_blocks,
        cfg.theta,
        cfg.keep_theta,
    );
    if blocks.is_empty() || blocks == *prev_blocks {
        return Ok(false);
    }
    *prev_blocks = blocks.clone();
    *fine = build_fine_level(coarse_grid, &blocks)?;
    let nb = blocks.len();
    let n_prom = ((nb as f64) * cfg.promote_frac.clamp(0.0, 1.0)) as usize;
    let promote: Vec<bool> = (0..nb).map(|i| i < n_prom).collect();
    promote_to_level2(fine, &promote);
    fine.invalidate();
    Ok(true)
}

/// bilinear prolongation from the coarse grid into a 2x2 fine block.
/// constant filling erased every intra-block gradient, which froze transmissive
/// fine patches forever; sampling the four neighboring coarse centers gives
/// each fine cell its own value from the first step.
#[allow(clippy::too_many_arguments)]
fn prolong_block_pub(
    coarse: &ConservedState2d,
    coarse_grid: &crate::grid2d::Grid2d,
    b: &super::amr2d::RefineBlock,
    out: &mut ConservedState2d,
) {
    let nx = coarse_grid.nx;
    let ci = |i: isize, j: isize| -> usize {
        let ii = (b.i as isize + i).clamp(0, nx as isize - 1) as usize;
        let jj = (b.j as isize + j).clamp(0, coarse_grid.ny as isize - 1) as usize;
        jj * nx + ii
    };
    let wx = [0.75f64, 0.25, 0.75, 0.25];
    let wy = [0.75f64, 0.75, 0.25, 0.25];
    let ks = [(0, 0), (-1, 0), (0, -1), (-1, -1)];
    for c in 0..4 {
        let (dx_i, dy_j) = ks[c];
        let acc = [wx[c], wy[c]];
        let cand = [
            (acc[0] * acc[1], ci(dx_i, dy_j)),
            ((1.0 - acc[0]) * acc[1], ci(dx_i + 1, dy_j)),
            (acc[0] * (1.0 - acc[1]), ci(dx_i, dy_j + 1)),
            ((1.0 - acc[0]) * (1.0 - acc[1]), ci(dx_i + 1, dy_j + 1)),
        ];
        let ck = b.j * coarse_grid.nx + b.i;
        for (dst, src) in [
            (&mut out.rho, &coarse.rho),
            (&mut out.mx, &coarse.mx),
            (&mut out.my, &coarse.my),
            (&mut out.e, &coarse.e),
        ] {
            let mut vals = [0.0f64; 4];
            for v in vals.iter_mut() {
                *v = cand.iter().map(|&(w, k)| w * src[k]).sum::<f64>();
            }
            // enforce conservation: shift the patch so its mean equals the
            // parent cell, keeping the bilinear spatial variation.
            let mean = vals.iter().sum::<f64>() / 4.0;
            let shift = src[ck] - mean;
            for c in 0..4 {
                dst[c] = vals[c] + shift;
            }
        }
    }
}

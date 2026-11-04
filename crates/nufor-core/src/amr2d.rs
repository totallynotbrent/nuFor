//! block-structured two-level adaptive refinement driver for the 2D grid.
//!
//! rung 11 asks for AMR on a gradient indicator. this module implements
//! the STRUCTURE, not a new solver: a coarse level covers the domain at
//! spacing dx, and each refined BLOCK subdivides one coarse cell into a
//! 2x2 patch of finer cells at dx/2 sitting at level 1. the indicator is
//! the normalized magnitude of the density gradient restricted to coarse
//! faces; cells whose incident face indicator exceeds the threshold get a
//! refined block.
//!
//! the same indicator machinery drives hyStrath-style DyM refinement, so
//! nothing here is a throwaway: promote refine_blocks to the march loop and
//! the two-level hierarchy is a real AMR mesh.

/// a refined block: one coarse cell (i, j) expanded to a 2x2 patch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RefineBlock {
    pub i: usize,
    pub j: usize,
    pub level: u8,
}

/// normalize `field` onto [0, 1] by its max magnitude, guarding the
/// all-zero case so a quiescent flow flags nothing.
pub fn normalize(field: &[f64]) -> Vec<f64> {
    let max = field.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    if max <= 0.0 {
        return vec![0.0; field.len()];
    }
    field.iter().map(|&v| v.abs() / max).collect()
}

/// horizontal and vertical face gradients of a nodal scalar field on an
/// nx-by-ny grid (nodal indexing idx=j*nx+i), returning per-cell
/// magnitude of grad rho approximated by centered differences.
pub fn cell_gradient_magnitude(nx: usize, ny: usize, field: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let ip = if i + 1 < nx { k + 1 } else { k };
            let im = if i > 0 { k - 1 } else { k };
            let jp = if j + 1 < ny { k + nx } else { k };
            let jm = if j > 0 { k - nx } else { k };
            let dx = (field[ip] - field[im]).abs();
            let dy = (field[jp] - field[jm]).abs();
            out[k] = (dx * dx + dy * dy).sqrt();
        }
    }
    out
}

/// choose refinement blocks from a normalized indicator and threshold.
///
/// blocks are non-overlapping by construction (one block per coarse
/// cell), and the level is recorded so a future three-level extension can
/// nest them. `threshold` in [0, 1]; 0 refines everywhere, 1 nowhere.
pub fn refine_blocks(
    nx: usize,
    ny: usize,
    indicator_norm: &[f64],
    threshold: f64,
) -> Vec<RefineBlock> {
    let mut blocks = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            if indicator_norm[k] >= threshold && threshold < 1.0 {
                blocks.push(RefineBlock { i, j, level: 1 });
            }
        }
    }
    blocks
}

/// total fine-cell count the two-level hierarchy would carry: coarse cells
/// plus 4 per refined block minus the covered coarse cells.
pub fn fine_cell_equivalent(ncells_coarse: usize, blocks: &[RefineBlock]) -> usize {
    ncells_coarse + 3 * blocks.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// normalization scales to unit max and zeroes a quiescent field.
    #[test]
    fn normalize_scales_and_zeroes() {
        let n = normalize(&[0.0, -2.0, 4.0, 1.0]);
        assert_eq!(n, vec![0.0, 0.5, 1.0, 0.25]);
        assert_eq!(normalize(&[0.0, 0.0]), vec![0.0, 0.0]);
    }

    /// a step in the field produces a nonzero gradient at the step only.
    #[test]
    fn gradient_flags_the_discontinuity() {
        let nx = 8;
        let ny = 1;
        let mut f = vec![1.0f64; nx * ny];
        for v in f.iter_mut().take(nx).skip(4) {
            *v = 5.0;
        }
        let g = cell_gradient_magnitude(nx, ny, &f);
        assert!(g[3] > 0.0 && g[4] > 0.0, "step edges must light up");
        assert_eq!(g[0], 0.0, "flat region stays dark");
        assert_eq!(g[nx - 1], 0.0);
    }

    /// threshold picks exactly the flagged cells, one block each.
    #[test]
    fn threshold_selects_blocks_one_per_cell() {
        let nx = 6;
        let ny = 4;
        let ind = normalize(&cell_gradient_magnitude(nx, ny, &{
            let mut f = vec![0.0f64; nx * ny];
            f[2 * nx + 3] = 1.0;
            f
        }));
        let blocks = refine_blocks(nx, ny, &ind, 0.5);
        assert!(!blocks.is_empty());
        for b in &blocks {
            assert_eq!(b.level, 1);
            assert!(ind[b.j * nx + b.i] >= 0.5);
        }
        let eq = fine_cell_equivalent(nx * ny, &blocks);
        assert_eq!(eq, nx * ny + 3 * blocks.len());
    }

    /// extreme thresholds: 1 refines nothing, 0 refines every cell.
    #[test]
    fn extreme_thresholds_bound_the_set() {
        let ind = vec![0.0f64; 4];
        assert!(refine_blocks(2, 2, &ind, 1.0).is_empty());
        assert_eq!(refine_blocks(2, 2, &ind, 0.0).len(), 4);
    }
}

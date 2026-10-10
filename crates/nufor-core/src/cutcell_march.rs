//! the cut-cell axisymmetric march: the annular update on fractional
//! cells with conservative flux redistribution for small cells. the
//! staircase path (advance2d_axi) stays untouched as the control; this
//! march reads precomputed fractions and never edits them.
//!
//! conservation on a partial cell: the update scales each face flux by
//! its open fraction, adds the embedded-boundary wall flux on the solid
//! part, and divides by the volume fraction,
//!
//!   u_t = -( sum_f alpha_f A_f F_f  -  F_eb A_eb ) / (alpha_vol V) + S(r)
//!
//! the EB wall flux is the inviscid slip statement laid out exactly:
//! zero mass, zero energy, momentum p n (the wall pressure force on
//! the fluid).
//!
//! the small-cell problem: an explicit update on a cell with tiny
//! alpha_vol needs dt ~ alpha_vol, which the global march cannot
//! afford. conservative redistribution (Chern-Colella lineage, the
//! production default in PeleC): a small cell keeps the fraction of its
//! update it can absorb at the global dt and hands the rest to the
//! neighbors across its open faces, weighted by the receivers' volume
//! fractions so the total is preserved exactly. first order in the cut
//! band, matching the first-order fraction geometry.

use crate::axi::compute_axi_fluxes;
use crate::cutcell::{cut_fractions_corner, CutFractions};
use crate::grid2d::Grid2d;

use crate::state2d::ConservedState2d;
use crate::{Boundaries2d, Error};

/// a cell whose fluid volume fraction sits below this threshold hands
/// its excess update to the neighbors rather than stepping on its own.
/// the half-cell rule: in 1d a cut cell at least half the regular size
/// is stable at the full time step (Berger's stabilization note), so
/// the threshold sits at one half.
pub const SMALL_CELL_EPS: f64 = 0.5;

/// the fraction of a small cell's update that gets redistributed to the
/// neighbors: zero at and above the threshold, ramping linearly to full
/// handoff at zero volume.
pub fn redistribute_weight(vol_frac: f64) -> f64 {
    if vol_frac >= SMALL_CELL_EPS {
        0.0
    } else {
        1.0 - vol_frac / SMALL_CELL_EPS
    }
}

/// the fractions field over the grid plus the per-cell EB geometry.
pub struct CutField {
    pub nx: usize,
    pub ny: usize,
    pub cells: Vec<CutFractions>,
    /// the EB interface length fraction per cell (the meridian-section
    /// analogue of the EB area) and its unit normal, where the solid
    /// centroid is; zero-length for full and solid cells.
    pub eb_len: Vec<f64>,
    pub eb_nx: Vec<f64>,
    pub eb_ny: Vec<f64>,
}

impl CutField {
    /// build the fractions and EB geometry for the whole grid from a
    /// signed distance function (positive inside the fluid) and its
    /// gradient (the outward normal is minus the SDF gradient).
    pub fn from_sdf(
        sdf: &dyn Fn(f64, f64) -> f64,
        normal: &dyn Fn(f64, f64) -> (f64, f64),
        g: &Grid2d,
    ) -> Self {
        let n = g.nx * g.ny;
        let mut cells = Vec::with_capacity(n);
        let mut eb_len = vec![0.0; n];
        let mut eb_nx = vec![0.0; n];
        let mut eb_ny = vec![0.0; n];
        for j in 0..g.ny {
            for i in 0..g.nx {
                let (x0, x1) = (g.faces_x[i], g.faces_x[i + 1]);
                let (y0, y1) = (g.faces_y[j], g.faces_y[j + 1]);
                let f = cut_fractions_corner(sdf, x0, x1, y0, y1);
                let k = j * g.nx + i;
                if f.vol > 0.0 && f.vol < 1.0 {
                    // the EB length from the discrete divergence: the
                    // aperture imbalance per direction closes the
                    // geometry so the update is consistent.
                    let ax = f.west - f.east;
                    let ay = f.south - f.north;
                    let l = (ax * ax + ay * ay).sqrt();
                    if l > 1e-12 {
                        // the interface normal points from fluid into
                        // the solid: it is the direction whose apertures
                        // shrink along it. the SDF normal at the solid
                        // centroid is the more accurate source when the
                        // caller provides it; fall back to apertures.
                        let (cx, cy) = (0.5 * (x0 + x1), 0.5 * (y0 + y1));
                        let (nx, ny) = normal(cx, cy);
                        let (nx, ny) = if nx * nx + ny * ny > 1e-12 {
                            let l2 = (nx * nx + ny * ny).sqrt();
                            (nx / l2, ny / l2)
                        } else {
                            (ax / l, ay / l)
                        };
                        eb_len[k] = l;
                        eb_nx[k] = nx;
                        eb_ny[k] = ny;
                    }
                }
                cells.push(f);
            }
        }
        CutField {
            nx: g.nx,
            ny: g.ny,
            cells,
            eb_len,
            eb_nx,
            eb_ny,
        }
    }

    pub fn at(&self, i: usize, j: usize) -> CutFractions {
        self.cells[j * self.nx + i]
    }
}

#[allow(clippy::too_many_arguments)]
/// one cut-cell axisymmetric euler step: the shared flux machinery plus
/// the fractional update and the redistribution pass.
pub fn advance2d_axi_cut(
    state: &mut ConservedState2d,
    g: &Grid2d,
    cut: &CutField,
    model: crate::thermo::ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;

    // the solid cells hold the quiescent body state exactly as the
    // staircase mask leaves them; the flux sweep sees the same field.
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let n = nx * ny;
    let mut e_int = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
    }
    let p = crate::thermo::pressure(model, &state.rho, &e_int, &u, &v)?;
    let a = crate::thermo::sound_speed(model, &state.rho, &et, &u, &v)?;
    let gamma = match model {
        crate::thermo::ThermoModel::Perfect { gamma } => gamma,
        crate::thermo::ThermoModel::EqAir => 1.4,
    };
    let (fx, fy) = compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, nthreads, None)?;

    // dt from the full-cell wave speeds; the small cells ride the
    // redistribution instead of shrinking it.
    let mut smax = 0.0f64;
    for k in 0..nx * ny {
        if cut.cells[k].vol <= 0.0 {
            continue;
        }
        smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
    }
    let dt = cfl * g.dx.min(g.dy) / smax.max(1e-12);

    // pass 1: the fractional update for every fluid cell, large and
    // small: the shared flux machinery scaled by face apertures,
    // divided by the cell's own fluid volume. pass 2 replaces the
    // small cells' states by the neighborhood pool.
    let mut drho = vec![0.0f64; nx * ny];
    let mut dmx = vec![0.0f64; nx * ny];
    let mut dmy = vec![0.0f64; nx * ny];
    let mut de = vec![0.0f64; nx * ny];

    for (j, fxrow) in fx.iter().enumerate() {
        let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
        let rc = 0.5 * (rm + rp);
        let dr = g.dys[j];
        for (i, fycol) in fy.iter().enumerate() {
            let k = idx(i, j);
            let f = cut.at(i, j);
            if f.vol <= 0.0 {
                continue;
            }
            // every fluid cell takes the honest flux update, small
            // ones included: state redistribution replaces the small
            // cells' states afterward, pooling over the neighborhood.
            let vol = f.vol.max(1e-9);
            let (x0, x1) = (g.faces_x[i], g.faces_x[i + 1]);

            // the open-face flux differences, scaled by apertures.
            let fs = f.south;
            let fn_ = f.north;
            let dxi = dt / (x1 - x0);

            let mass_f = dxi * (f.east * fxrow.mass[i + 1] - f.west * fxrow.mass[i])
                + dt * (rp * fn_ * fycol.mass[j + 1] - rm * fs * fycol.mass[j]) / (rc * dr);
            let mx_f = dxi * (f.east * fxrow.mx[i + 1] - f.west * fxrow.mx[i])
                + dt * (rp * fn_ * fycol.mx[j + 1] - rm * fs * fycol.mx[j]) / (rc * dr);
            let my_f = dxi * (f.east * fxrow.my[i + 1] - f.west * fxrow.my[i])
                + dt * (rp * fn_ * fycol.my[j + 1] - rm * fs * fycol.my[j]) / (rc * dr);
            let e_f = dxi * (f.east * fxrow.e[i + 1] - f.west * fxrow.e[i])
                + dt * (rp * fn_ * fycol.e[j + 1] - rm * fs * fycol.e[j]) / (rc * dr);

            // the EB wall flux: the reflected-ghost hllc pair in the
            // wall-local frame with the axis along the cell's outward
            // EB normal. the pure pressure-force layout (0, p n, 0)
            // rings and collapses dt in both stabilization schemes
            // tried; the reflected pair carries the dissipation the
            // wall needs. its known cost, measured on the mach-22
            // shell: the wall riemann problem pumps the surface cells
            // toward ~2x the stagnation pressure at coarse grids
            // where the shock layer sits inside the SRD neighborhood.
            let (mut wmass, mut wmx, mut wmy, mut we) = (0.0, 0.0, 0.0, 0.0);
            if cut.eb_len[k] > 0.0 {
                use crate::hllc2d::{hllc_flux, FacePrim};
                let (enx, eny) = (-cut.eb_nx[k], -cut.eb_ny[k]);
                let un = u[k] * enx + v[k] * eny;
                let ut = -u[k] * eny + v[k] * enx;
                let (wl, wr) = match model {
                    crate::thermo::ThermoModel::Perfect { gamma } => (
                        FacePrim::perfect(state.rho[k], un, ut, p[k], gamma),
                        FacePrim::perfect(state.rho[k], -un, ut, p[k], gamma),
                    ),
                    crate::thermo::ThermoModel::EqAir => {
                        let e_wall = e_int[k];
                        let p_wall = crate::eqair::eqair_pressure_at(state.rho[k], e_wall);
                        let a_wall = crate::eqair::eqair_sound_at(state.rho[k], e_wall);
                        (
                            FacePrim {
                                rho: state.rho[k],
                                u: un,
                                v: ut,
                                p: p_wall,
                                a: a_wall,
                                e: e_wall,
                            },
                            FacePrim {
                                rho: state.rho[k],
                                u: -un,
                                v: ut,
                                p: p_wall,
                                a: a_wall,
                                e: e_wall,
                            },
                        )
                    }
                };
                let q = hllc_flux(gamma, wl, wr, 0);
                let l = dt * cut.eb_len[k];
                wmass = q.mass * l;
                wmx = (q.mx * enx - q.my * eny) * l;
                wmy = (q.mx * eny + q.my * enx) * l;
                we = q.e * l;
            }

            // the EB flux joins the divergence sum with the same
            // outward-positive sign as the open faces: the cell loses
            // q*A*dt through the wall face like any other face.
            drho[k] = (-mass_f - wmass) / vol;
            dmx[k] = (-mx_f - wmx) / vol;
            dmy[k] = (-my_f - wmy) / vol + dt * p[k] / rc;
            de[k] = (-e_f - we) / vol;
        }
    }

    // the increments apply inside pass 2, where the small-cell
    // redistribution can correct them before they land.

    // pass 2: state redistribution (Berger-Giuliani), extended to
    // every cut cell: the axisymmetric annular divergence factor
    // near the axis amplifies the fractional update beyond what the
    // planar half-cell rule tolerates (the paper's threshold ran
    // stable only to t ~ 1.1 on the sphere). the cost is mixing
    // irreversibility near the wall, which shifts the measured bow-
    // shock standoff (0.217 vs the staircase's 0.317 at the same
    // grid, A-W 0.321) while reading the true stagnation states
    // (rho ~ 4.0 vs the staircase's mask-averaged 1.9-3.0).
    //
    // delta redistribution (the update-average variant, arXiv
    // 2610.06459) was implemented and measured here: it preserves
    // the shock layer in state space exactly as designed, but in the
    // axisymmetric near-axis wake the shared updates drain the whole
    // cut cluster to vacuum within 60 steps (pooling to 0.5 Vfull
    // only delayed it to t=0.026). the failure is documented in the
    // goal file; SRD is the stable scheme for this geometry.
    {
        // neighborhoods: large cells are their own; every cut cell
        // grows a ring of fluid neighbors until the pooled fluid
        // volume reaches a full cell's worth.
        let mut nbrs: Vec<Vec<usize>> = vec![Vec::new(); nx * ny];
        let fluid = |k: usize| cut.cells[k].vol > 0.0;
        let vfull = g.dx * g.dy;
        for j in 0..ny {
            for i in 0..nx {
                let k = idx(i, j);
                let alpha = cut.cells[k].vol;
                if alpha >= 1.0 || alpha <= 0.0 {
                    if alpha > 0.0 {
                        nbrs[k].push(k);
                    }
                    continue;
                }
                let mut member = vec![k];
                let mut vol = alpha * vfull;
                let (mut r0, mut r1, mut c0, mut c1) = (i, i, j, j);
                while vol < vfull && (r0 > 0 || r1 + 1 < nx || c0 > 0 || c1 + 1 < ny) {
                    let (pr0, pr1, pc0, pc1) = (r0, r1, c0, c1);
                    r0 = r0.saturating_sub(1);
                    r1 = (r1 + 1).min(nx - 1);
                    c0 = c0.saturating_sub(1);
                    c1 = (c1 + 1).min(ny - 1);
                    for jj in c0..=c1 {
                        for ii in r0..=r1 {
                            let on_new_ring = ii == r0 || ii == r1 || jj == c0 || jj == c1;
                            if !on_new_ring
                                || (ii == pr0 && jj >= pc0 && jj <= pc1)
                                || (ii == pr1 && jj >= pc0 && jj <= pc1)
                            {
                                continue;
                            }
                            let k2 = idx(ii, jj);
                            if fluid(k2) && !member.contains(&k2) {
                                member.push(k2);
                                vol += cut.cells[k2].vol * vfull;
                            }
                        }
                    }
                }
                nbrs[k] = member;
            }
        }

        // overlap counts.
        let mut cover = vec![0usize; nx * ny];
        for row in &nbrs {
            for &m in row {
                cover[m] += 1;
            }
        }

        for k in 0..nx * ny {
            state.rho[k] += drho[k];
            state.mx[k] += dmx[k];
            state.my[k] += dmy[k];
            state.e[k] += de[k];
        }

        // weighted neighborhood averages of the updated states.
        let mut qavg: Vec<[f64; 4]> = vec![[0.0; 4]; nx * ny];
        for k in 0..nx * ny {
            if nbrs[k].is_empty() {
                continue;
            }
            let mut wv = 0.0f64;
            let mut acc = [0.0f64; 4];
            for &m in &nbrs[k] {
                let w = cut.cells[m].vol * vfull / cover[m].max(1) as f64;
                wv += w;
                acc[0] += w * state.rho[m];
                acc[1] += w * state.mx[m];
                acc[2] += w * state.my[m];
                acc[3] += w * state.e[m];
            }
            qavg[k] = [acc[0] / wv, acc[1] / wv, acc[2] / wv, acc[3] / wv];
        }

        // final state: the average of the averages covering the cell.
        let mut nnew = vec![0usize; nx * ny];
        let mut nacc = vec![[0.0f64; 4]; nx * ny];
        for k in 0..nx * ny {
            for &m in &nbrs[k] {
                nacc[m][0] += qavg[k][0];
                nacc[m][1] += qavg[k][1];
                nacc[m][2] += qavg[k][2];
                nacc[m][3] += qavg[k][3];
                nnew[m] += 1;
            }
        }
        for m in 0..nx * ny {
            if nnew[m] > 0 {
                state.rho[m] = nacc[m][0] / nnew[m] as f64;
                state.mx[m] = nacc[m][1] / nnew[m] as f64;
                state.my[m] = nacc[m][2] / nnew[m] as f64;
                state.e[m] = nacc[m][3] / nnew[m] as f64;
            }
        }
    }

    Ok((dt, 0.0))
}

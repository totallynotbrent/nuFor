//! spalart-allmaras transport step: upwind advection, central diffusion, and
//! the production/destruction source, advancing nu_tilde over one dt.
//!
//! advection is first-order upwind on the cell-centered field using the
//! face-averaged mean-flow velocity (the standard treatment for turbulence
//! scalars); diffusion is the central two-point flux with the sa diffusivity
//! (nu + nu_tilde)/sigma plus the c_b2/sigma cross term; the source is the
//! 1994 production/destruction pair. the update advances the conservative
//! rho*nu_tilde and divides by the cell density, so a uniform field in
//! uniform flow is exactly preserved. solid walls zero nu_tilde in the ghost
//! layer (the sa wall condition); open sides copy the interior.

use crate::error::Error;
use crate::grid2d::Grid2d;
use crate::sa::{source, C_B2, SIGMA};
use crate::solver2d::Boundaries2d;
use crate::state2d::{cons_to_prim2d, ConservedState2d};

/// the sa model state: the transported field, the wall distance it needs,
/// and the model parameters (freestream value, prandtl numbers, and the
/// molecular viscosity the closures use).
#[derive(Debug, Clone, PartialEq)]
pub struct TurbState {
    /// the modified turbulent viscosity nu_tilde per cell.
    pub nu_tilde: Vec<f64>,
    /// per-cell wall distance feeding the destruction term.
    pub d: Vec<f64>,
    /// the model parameters.
    pub params: SaParams,
}

/// the sa model parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SaParams {
    /// molecular dynamic viscosity.
    pub mu: f64,
    /// molecular prandtl number.
    pub pr: f64,
    /// freestream nu_tilde, the far-field/inflow value (convention 3*nu).
    pub nu_tilde_inf: f64,
    /// turbulent prandtl number for the eddy heat flux.
    pub pr_t: f64,
}

/// vorticity magnitude s = sqrt(2 omega_ij omega_ij) from centered velocity
/// gradients; in 2d this is |dv/dx - du/dy|. on stretched axes the
/// derivatives use the true neighbor distances (the plain central
/// difference when uniform).
fn vorticity(u: &[f64], v: &[f64], g: &Grid2d, unif_x: bool, unif_y: bool) -> Vec<f64> {
    let (dx, dy) = (g.dx, g.dy);
    let (nx, ny) = (g.nx, g.ny);
    let w = nx + 2;
    let pad = |f: &[f64]| -> Vec<f64> {
        let mut out = vec![0.0; w * (ny + 2)];
        for j in 0..ny {
            for i in 0..nx {
                out[(j + 1) * w + (i + 1)] = f[j * nx + i];
            }
        }
        for j in 0..ny {
            out[(j + 1) * w] = out[(j + 1) * w + 1];
            out[(j + 1) * w + nx + 1] = out[(j + 1) * w + nx];
        }
        for i in 0..w {
            out[i] = out[i + w];
            out[(ny + 1) * w + i] = out[ny * w + i];
        }
        out
    };
    let pu = pad(u);
    let pv = pad(v);
    let mut out = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let (ip, jp) = (i + 1, j + 1);
            let c = j * nx + i;
            let dvdx = if unif_x {
                (pv[jp * w + (i + 2)] - pv[jp * w + i]) / (2.0 * dx)
            } else {
                // the 3-point lagrange derivative at the cell center over
                // its neighbors, on their true positions.
                let (xl, xr) = if i == 0 {
                    (2.0 * g.faces_x[0] - g.faces_x[1], g.faces_x[1])
                } else if i == nx - 1 {
                    (g.faces_x[nx - 1], 2.0 * g.faces_x[nx] - g.faces_x[nx - 1])
                } else {
                    (g.faces_x[i], g.faces_x[i + 1])
                };
                let xc = 0.5 * (g.faces_x[i] + g.faces_x[i + 1]);
                let (vl, vc, vr) = (pv[jp * w + i], pv[jp * w + ip], pv[jp * w + (i + 2)]);
                let w_l = (xc - xr) / ((xl - xc) * (xl - xr));
                let w_c = (2.0 * xc - xl - xr) / ((xc - xl) * (xc - xr));
                let w_r = (xc - xl) / ((xr - xl) * (xr - xc));
                vl * w_l + vc * w_c + vr * w_r
            };
            let dudy = if unif_y {
                (pu[(j + 2) * w + ip] - pu[j * w + ip]) / (2.0 * dy)
            } else {
                let (yl, yu) = if j == 0 {
                    (2.0 * g.faces_y[0] - g.faces_y[1], g.faces_y[1])
                } else if j == ny - 1 {
                    (g.faces_y[ny - 1], 2.0 * g.faces_y[ny] - g.faces_y[ny - 1])
                } else {
                    (g.faces_y[j], g.faces_y[j + 1])
                };
                let yc = 0.5 * (g.faces_y[j] + g.faces_y[j + 1]);
                let (vl, vc, vu) = (pu[j * w + ip], pu[jp * w + ip], pu[(j + 2) * w + ip]);
                let w_l = (yc - yu) / ((yl - yc) * (yl - yu));
                let w_c = (2.0 * yc - yl - yu) / ((yc - yl) * (yc - yu));
                let w_u = (yc - yl) / ((yu - yl) * (yu - yc));
                vl * w_l + vc * w_c + vu * w_u
            };
            out[c] = (dvdx - dudy).abs();
        }
    }
    out
}

/// true when a boundary side is a solid wall (sa zeroes nu_tilde there).
fn is_wall(bc: &crate::solver2d::Bc2d) -> bool {
    matches!(
        bc,
        crate::solver2d::Bc2d::SlipWall | crate::solver2d::Bc2d::NoSlipWall
    )
}

/// padded nu_tilde with wall-aware ghosts: solid sides zero, open sides copy.
fn pad_scalar(nt: &[f64], nx: usize, ny: usize, bc: &Boundaries2d) -> Vec<f64> {
    let w = nx + 2;
    let mut out = vec![0.0; w * (ny + 2)];
    for j in 0..ny {
        for i in 0..nx {
            out[(j + 1) * w + (i + 1)] = nt[j * nx + i];
        }
    }
    for j in 0..ny {
        out[(j + 1) * w] = if is_wall(&bc.west) {
            0.0
        } else {
            out[(j + 1) * w + 1]
        };
        out[(j + 1) * w + nx + 1] = if is_wall(&bc.east) {
            0.0
        } else {
            out[(j + 1) * w + nx]
        };
    }
    for i in 0..w {
        out[i] = if is_wall(&bc.south) { 0.0 } else { out[i + w] };
        out[(ny + 1) * w + i] = if is_wall(&bc.north) {
            0.0
        } else {
            out[ny * w + i]
        };
    }
    out
}

/// advance nu_tilde by dt: upwind advection + central diffusion + source.
///
/// the wall distance and molecular viscosity come from the turb state.
pub fn advance_turb(
    turb: &mut TurbState,
    state: &ConservedState2d,
    g: &Grid2d,
    bc: &Boundaries2d,
    dt: f64,
) -> Result<(), Error> {
    let mu_lam = turb.params.mu;
    let d = &turb.d;
    let nx = g.nx;
    let ny = g.ny;
    let n = nx * ny;
    if turb.nu_tilde.len() != n || d.len() != n {
        return Err(Error::InvalidArgs);
    }
    if dt <= 0.0 {
        return Ok(());
    }
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    // per-cell width ratios; uniform axes keep the scalar fast path.
    let unif_x = (0..nx).all(|i| (g.dxs[i] - g.dxs[0]).abs() <= 1e-9 * g.dxs[0].abs().max(1e-12));
    let unif_y = (0..ny).all(|j| (g.dys[j] - g.dys[0]).abs() <= 1e-9 * g.dys[0].abs().max(1e-12));
    let s = vorticity(&u, &v, g, unif_x, unif_y);
    let pn = pad_scalar(&turb.nu_tilde, nx, ny, bc);
    let w = nx + 2;
    let mut new_nt = turb.nu_tilde.clone();
    let inv_xs: Vec<f64> = g.dxs.iter().map(|w| 1.0 / w).collect();
    let inv_ys: Vec<f64> = g.dys.iter().map(|w| 1.0 / w).collect();
    let (inv_x, inv_y) = (
        if unif_x { 1.0 / g.dx } else { 0.0 },
        if unif_y { 1.0 / g.dy } else { 0.0 },
    );

    for j in 0..ny {
        for i in 0..nx {
            let c = j * nx + i;
            let rho = state.rho[c].max(1e-12);
            let nu_c = mu_lam / rho;
            // face fluxes on the four sides of the cell; net update is the
            // difference (outflow minus inflow) for advection, and (in minus
            // out) for the laplacian part of the diffusion.
            let mut f_adv = [0.0f64; 4];
            let mut f_dif = [0.0f64; 4];
            // x-faces: face f sits between padded columns f and f+1; cell i owns
            // faces i (left) and i+1 (right). face values average the two
            // adjacent cells (single-sided at the domain edges).
            for (k, f) in [i, i + 1].iter().enumerate() {
                let f = *f;
                let (la, lb) = (f, f + 1);
                let (cell_l, cell_r) = (f as isize - 1, f as isize);
                let u_f = if f == 0 {
                    u[j * nx]
                } else if f == nx {
                    u[j * nx + nx - 1]
                } else {
                    0.5 * (u[j * nx + (cell_l as usize)] + u[j * nx + (cell_r as usize)])
                };
                let rho_f = if f == 0 {
                    state.rho[j * nx]
                } else if f == nx {
                    state.rho[j * nx + nx - 1]
                } else {
                    0.5 * (state.rho[j * nx + (cell_l as usize)]
                        + state.rho[j * nx + (cell_r as usize)])
                };
                let left = pn[(j + 1) * w + la];
                let right = pn[(j + 1) * w + lb];
                let nt_up = if u_f >= 0.0 { left } else { right };
                f_adv[k] = rho_f * u_f * nt_up;
                let nt_f = 0.5 * (left + right);
                let nu_f = mu_lam / rho_f.max(1e-12);
                let inv = if unif_x {
                    inv_x
                } else {
                    // the true center-to-center distance across face f.
                    let dc = if f == 0 {
                        g.centers_x[j * nx] - g.faces_x[0]
                    } else if f == nx {
                        g.faces_x[nx] - g.centers_x[j * nx + nx - 1]
                    } else {
                        g.centers_x[j * nx + f] - g.centers_x[j * nx + f - 1]
                    };
                    1.0 / dc
                };
                f_dif[k] = (nu_f + nt_f) / SIGMA * (right - left) * inv;
            }
            // y-faces: face f sits between padded rows f and f+1; cell j owns
            // faces j (bottom) and j+1 (top).
            for (k, f) in [j, j + 1].iter().enumerate() {
                let f = *f;
                let (ra, rb) = (f, f + 1);
                let (cell_b, cell_t) = (f as isize - 1, f as isize);
                let v_f = if f == 0 {
                    v[i]
                } else if f == ny {
                    v[(ny - 1) * nx + i]
                } else {
                    0.5 * (v[(cell_b as usize) * nx + i] + v[(cell_t as usize) * nx + i])
                };
                let rho_f = if f == 0 {
                    state.rho[i]
                } else if f == ny {
                    state.rho[(ny - 1) * nx + i]
                } else {
                    0.5 * (state.rho[(cell_b as usize) * nx + i]
                        + state.rho[(cell_t as usize) * nx + i])
                };
                let below = pn[ra * w + i + 1];
                let above = pn[rb * w + i + 1];
                let nt_up = if v_f >= 0.0 { below } else { above };
                f_adv[k + 2] = rho_f * v_f * nt_up;
                let nt_f = 0.5 * (below + above);
                let nu_f = mu_lam / rho_f.max(1e-12);
                let inv = if unif_y {
                    inv_y
                } else {
                    // the true center-to-center distance across face f.
                    let dc = if f == 0 {
                        g.centers_y[i] - g.faces_y[0]
                    } else if f == ny {
                        g.faces_y[ny] - g.centers_y[(ny - 1) * nx + i]
                    } else {
                        g.centers_y[f * nx + i] - g.centers_y[(f - 1) * nx + i]
                    };
                    1.0 / dc
                };
                f_dif[k + 2] = (nu_f + nt_f) / SIGMA * (above - below) * inv;
            }
            // conservative flux differences: advection out minus in, diffusion
            // in minus out (the laplacian has a plus sign in the sa equation).
            let inv_x_i = if unif_x { inv_x } else { inv_xs[i] };
            let inv_y_j = if unif_y { inv_y } else { inv_ys[j] };
            let adv_net = (f_adv[1] - f_adv[0]) * inv_x_i + (f_adv[3] - f_adv[2]) * inv_y_j;
            let dif_net = (f_dif[1] - f_dif[0]) * inv_x_i + (f_dif[3] - f_dif[2]) * inv_y_j;
            // the c_b2/sigma cross term from cell-centered gradients.
            let (ip, jp) = (i + 1, j + 1);
            let (gx, gy) = if unif_x && unif_y {
                (
                    (pn[jp * w + i + 2] - pn[jp * w + i]) / (2.0 * g.dx),
                    (pn[(j + 2) * w + ip] - pn[j * w + ip]) / (2.0 * g.dy),
                )
            } else {
                // lagrange derivative at the cell center on true positions
                // (ghost centers reflected across the domain faces).
                let xg = if i == 0 {
                    2.0 * g.faces_x[0] - g.centers_x[c]
                } else {
                    g.centers_x[c - 1]
                };
                let xgr = if i == nx - 1 {
                    2.0 * g.faces_x[nx] - g.centers_x[c]
                } else {
                    g.centers_x[c + 1]
                };
                let yg = if j == 0 {
                    2.0 * g.faces_y[0] - g.centers_y[c]
                } else {
                    g.centers_y[c - nx]
                };
                let ygt = if j == ny - 1 {
                    2.0 * g.faces_y[ny] - g.centers_y[c]
                } else {
                    g.centers_y[c + nx]
                };
                let xc = g.centers_x[c];
                let yc = g.centers_y[c];
                let lag = |v0: f64, v1: f64, v2: f64, a: f64, b: f64, cc: f64| {
                    v0 * (b - cc) / ((a - b) * (a - cc))
                        + v1 * (2.0 * b - a - cc) / ((b - a) * (b - cc))
                        + v2 * (b - a) / ((cc - a) * (cc - b))
                };
                (
                    lag(
                        pn[jp * w + i],
                        pn[jp * w + ip],
                        pn[jp * w + (i + 2)],
                        xg,
                        xc,
                        xgr,
                    ),
                    lag(
                        pn[j * w + ip],
                        pn[jp * w + ip],
                        pn[(j + 2) * w + ip],
                        yg,
                        yc,
                        ygt,
                    ),
                )
            };
            let cb2_term = C_B2 / SIGMA * (gx * gx + gy * gy);
            let src = source(turb.nu_tilde[c], s[c], d[c], nu_c, dt, rho);
            // advection is conservative in rho*nu_tilde (divide by rho);
            // diffusion is kinematic, acting on nu_tilde directly.
            let rhs = -adv_net / rho + dif_net + cb2_term + src;
            new_nt[c] = turb.nu_tilde[c] + dt * rhs;
        }
    }
    // positivity: the transported field is physically non-negative.
    for x in new_nt.iter_mut() {
        *x = x.max(0.0);
    }
    turb.nu_tilde = new_nt;
    Ok(())
}

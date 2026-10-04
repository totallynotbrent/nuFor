//! the 3d viscous operator: navier-stokes diffusive fluxes in three
//! directions on the uniform structured grid, with the no-slip wall
//! machinery the operator needs on any of the six faces.
//!
//! the layout matches solver3d: idx = (k*ny + j)*nx + i. interior gradients
//! are centered; wall-adjacent cells use the quadratic-consistent one-sided
//! stencils the 2d operator uses, so the wall flux and the wall-cell
//! balance agree exactly.

use crate::eos3d::eos_pressure3d;
use crate::grid3d::Grid3d;
use crate::solver3d::advance3d;
use crate::state3d::{cons_to_prim3d, ConservedState3d};

/// which of the six domain faces are solid; every open face stays
/// transmissive (matching the inviscid step's treatment).
#[derive(Debug, Clone, Copy, Default)]
pub struct Walls3d {
    pub west: bool,
    pub east: bool,
    pub south: bool,
    pub north: bool,
    pub bottom: bool,
    pub top: bool,
}

/// the physical parameters of the viscous fluid (3d mirror of ViscParams).
#[derive(Debug, Clone, Copy)]
pub struct ViscParams3d {
    pub mu: f64,
    pub pr: f64,
}

/// the (u, v, w, t) gradients at one cell.
#[derive(Debug, Clone, Copy, Default)]
struct Grads3d {
    ux: f64,
    uy: f64,
    uz: f64,
    vx: f64,
    vy: f64,
    vz: f64,
    wx: f64,
    wy: f64,
    wz: f64,
    tx: f64,
    ty: f64,
    tz: f64,
}

/// the one-sided wall gradient the 2d operator proved: exact for a profile
/// that is quadratic near the wall. used at solid faces in the axis normal
/// to the wall, for the wall-parallel velocity and temperature.
fn wall_grad(p0: f64, p1: f64, d: f64) -> f64 {
    (9.0 * p0 - p1) / (3.0 * d)
}

/// the wall-cell's own normal gradient so the cell update balances the wall
/// flux: the 2d wall_cell_grad form, (3 p0 + p1) / (3 d), signed by flip.
fn wall_cell_grad_val(p0: f64, p1: f64, d: f64, flip: f64) -> f64 {
    flip * (3.0 * p0 + p1) / (3.0 * d)
}

/// centered gradient when both neighbors exist, one-sided low at a low-side
/// wall (mirror value at the wall: the gradient through the ghost pair).
fn grad_axis(f: &[f64], lo: Option<usize>, hi: Option<usize>, c: usize, d: f64, flip: f64) -> f64 {
    match (lo, hi) {
        (Some(a), Some(b)) => (f[b] - f[a]) / (2.0 * d),
        // low-side wall: the ghost mirrors the center, so the centered
        // stencil degrades to the one-sided wall form.
        (None, Some(b)) => flip * (f[b] - f[c]) / d,
        (Some(a), None) => -flip * (f[c] - f[a]) / d,
        _ => 0.0,
    }
}

/// build the per-cell gradients, with the wall cells' normal components
/// replaced by the wall-consistent one-sided values.
#[allow(clippy::too_many_arguments)]
fn cell_grads3d(
    u: &[f64],
    v: &[f64],
    w: &[f64],
    t: &[f64],
    g: &Grid3d,
    walls: &Walls3d,
) -> Vec<Grads3d> {
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let n = nx * ny * nz;
    let idx = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
    let mut out = vec![Grads3d::default(); n];
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let c = idx(i, j, k);
                // neighbors per axis: solid walls remove the ghost side
                // (the one-sided forms then apply); open faces keep the
                // transmissive mirror (the ghost copies the edge cell).
                let lo_x = if i == 0 && !walls.west {
                    Some(c)
                } else if i > 0 {
                    Some(idx(i - 1, j, k))
                } else {
                    None
                };
                let hi_x = if i == nx - 1 && !walls.east {
                    Some(c)
                } else if i < nx - 1 {
                    Some(idx(i + 1, j, k))
                } else {
                    None
                };
                let lo_y = if j == 0 && !walls.south {
                    Some(c)
                } else if j > 0 {
                    Some(idx(i, j - 1, k))
                } else {
                    None
                };
                let hi_y = if j == ny - 1 && !walls.north {
                    Some(c)
                } else if j < ny - 1 {
                    Some(idx(i, j + 1, k))
                } else {
                    None
                };
                let lo_z = if k == 0 && !walls.bottom {
                    Some(c)
                } else if k > 0 {
                    Some(idx(i, j, k - 1))
                } else {
                    None
                };
                let hi_z = if k == nz - 1 && !walls.top {
                    Some(c)
                } else if k < nz - 1 {
                    Some(idx(i, j, k + 1))
                } else {
                    None
                };
                // at solid walls the one-sided flip makes (b - c) carry the
                // wall's derivative direction (wall below: +, wall above: -).
                let (fx, fy, fz) = (
                    if i == 0 && walls.west {
                        1.0
                    } else if i == nx - 1 && walls.east {
                        -1.0
                    } else {
                        1.0
                    },
                    if j == 0 && walls.south {
                        1.0
                    } else if j == ny - 1 && walls.north {
                        -1.0
                    } else {
                        1.0
                    },
                    if k == 0 && walls.bottom {
                        1.0
                    } else if k == nz - 1 && walls.top {
                        -1.0
                    } else {
                        1.0
                    },
                );
                out[c] = Grads3d {
                    ux: grad_axis(u, lo_x, hi_x, c, g.dx, fx),
                    uy: grad_axis(u, lo_y, hi_y, c, g.dy, fy),
                    uz: grad_axis(u, lo_z, hi_z, c, g.dz, fz),
                    vx: grad_axis(v, lo_x, hi_x, c, g.dx, fx),
                    vy: grad_axis(v, lo_y, hi_y, c, g.dy, fy),
                    vz: grad_axis(v, lo_z, hi_z, c, g.dz, fz),
                    wx: grad_axis(w, lo_x, hi_x, c, g.dx, fx),
                    wy: grad_axis(w, lo_y, hi_y, c, g.dy, fy),
                    wz: grad_axis(w, lo_z, hi_z, c, g.dz, fz),
                    tx: grad_axis(t, lo_x, hi_x, c, g.dx, fx),
                    ty: grad_axis(t, lo_y, hi_y, c, g.dy, fy),
                    tz: grad_axis(t, lo_z, hi_z, c, g.dz, fz),
                };
                // wall cells: the normal gradient of the wall-parallel
                // velocity and temperature must balance the wall flux, so
                // replace it with the wall-consistent one-sided value.
                let mut gr = out[c];
                if i == 0 && walls.west {
                    gr.ux = wall_cell_grad_val(u[c], u[idx(1, j, k)], g.dx, 1.0);
                    gr.vx = wall_cell_grad_val(v[c], v[idx(1, j, k)], g.dx, 1.0);
                    gr.wx = wall_cell_grad_val(w[c], w[idx(1, j, k)], g.dx, 1.0);
                    gr.tx = wall_cell_grad_val(t[c], t[idx(1, j, k)], g.dx, 1.0);
                }
                if i == nx - 1 && walls.east {
                    let b = idx(nx - 2, j, k);
                    gr.ux = wall_cell_grad_val(u[c], u[b], g.dx, -1.0);
                    gr.vx = wall_cell_grad_val(v[c], v[b], g.dx, -1.0);
                    gr.wx = wall_cell_grad_val(w[c], w[b], g.dx, -1.0);
                    gr.tx = wall_cell_grad_val(t[c], t[b], g.dx, -1.0);
                }
                if j == 0 && walls.south {
                    let b = idx(i, 1, k);
                    gr.uy = wall_cell_grad_val(u[c], u[b], g.dy, 1.0);
                    gr.vy = wall_cell_grad_val(v[c], v[b], g.dy, 1.0);
                    gr.wy = wall_cell_grad_val(w[c], w[b], g.dy, 1.0);
                    gr.ty = wall_cell_grad_val(t[c], t[b], g.dy, 1.0);
                }
                if j == ny - 1 && walls.north {
                    let b = idx(i, ny - 2, k);
                    gr.uy = wall_cell_grad_val(u[c], u[b], g.dy, -1.0);
                    gr.vy = wall_cell_grad_val(v[c], v[b], g.dy, -1.0);
                    gr.wy = wall_cell_grad_val(w[c], w[b], g.dy, -1.0);
                    gr.ty = wall_cell_grad_val(t[c], t[b], g.dy, -1.0);
                }
                if k == 0 && walls.bottom {
                    let b = idx(i, j, 1);
                    gr.uz = wall_cell_grad_val(u[c], u[b], g.dz, 1.0);
                    gr.vz = wall_cell_grad_val(v[c], v[b], g.dz, 1.0);
                    gr.wz = wall_cell_grad_val(w[c], w[b], g.dz, 1.0);
                    gr.tz = wall_cell_grad_val(t[c], t[b], g.dz, 1.0);
                }
                if k == nz - 1 && walls.top {
                    let b = idx(i, j, nz - 2);
                    gr.uz = wall_cell_grad_val(u[c], u[b], g.dz, -1.0);
                    gr.vz = wall_cell_grad_val(v[c], v[b], g.dz, -1.0);
                    gr.wz = wall_cell_grad_val(w[c], w[b], g.dz, -1.0);
                    gr.tz = wall_cell_grad_val(t[c], t[b], g.dz, -1.0);
                }
                out[c] = gr;
            }
        }
    }
    out
}

/// the viscous flux across one interior x-face: the stress row tau_x
/// dotted with the face velocity plus heat conduction. returns
/// (x-momentum, y-momentum, z-momentum, energy) per unit area.
#[allow(clippy::too_many_arguments)]
fn visc_flux_x(
    ga: &Grads3d,
    gb: &Grads3d,
    ua: f64,
    va: f64,
    wa: f64,
    ub: f64,
    vb: f64,
    wb: f64,
    mu: f64,
    kappa: f64,
) -> (f64, f64, f64, f64) {
    // face-averaged gradients.
    let (ux, vx, wx, tx) = (
        0.5 * (ga.ux + gb.ux),
        0.5 * (ga.vx + gb.vx),
        0.5 * (ga.wx + gb.wx),
        0.5 * (ga.tx + gb.tx),
    );
    let (uy, _wy) = (0.5 * (ga.uy + gb.uy), 0.5 * (ga.wy + gb.wy));
    let (uz, _vz) = (0.5 * (ga.uz + gb.uz), 0.5 * (ga.vz + gb.vz));
    let (vy, wz) = (0.5 * (ga.vy + gb.vy), 0.5 * (ga.wz + gb.wz));
    let (uf, vf, wf) = (0.5 * (ua + ub), 0.5 * (va + vb), 0.5 * (wa + wb));
    // the stress row tau_x with the stokes div-velocity term.
    let divu = ux + vy + wz;
    let txx = mu * (2.0 * ux - 2.0 / 3.0 * divu);
    let txy = mu * (uy + vx);
    let txz = mu * (uz + wx);
    let q = kappa * tx;
    let fe = uf * txx + vf * txy + wf * txz + q;
    (txx, txy, txz, fe)
}
/// one explicit viscous substep added to the state (heun's first stage).
#[allow(clippy::too_many_arguments)]
pub fn add_viscous3d(
    state: &mut ConservedState3d,
    g: &Grid3d,
    gamma: f64,
    v: ViscParams3d,
    walls: &Walls3d,
    dt: f64,
) -> Result<(), crate::Error> {
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let n = nx * ny * nz;
    if state.rho.len() != n
        || state.mx.len() != n
        || state.my.len() != n
        || state.mz.len() != n
        || state.e.len() != n
    {
        return Err(crate::Error::InvalidArgs);
    }
    if v.mu <= 0.0 || v.pr <= 0.0 || dt <= 0.0 {
        return Ok(());
    }
    let (u, vv, ww, et) = cons_to_prim3d(&state.rho, &state.mx, &state.my, &state.mz, &state.e)?;
    let p = eos_pressure3d(gamma, &state.rho, &et, &u, &vv, &ww)?;
    let t: Vec<f64> = p.iter().zip(&state.rho).map(|(pp, r)| pp / r).collect();
    let grads = cell_grads3d(&u, &vv, &ww, &t, g, walls);
    let kappa = v.mu * gamma / ((gamma - 1.0) * v.pr);
    let mu = v.mu;
    let idx = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;

    // interior x-face fluxes: faces 1..nx-1 for every (j, k).
    for k in 0..nz {
        for j in 0..ny {
            for i in 1..nx {
                let (a, b) = (idx(i - 1, j, k), idx(i, j, k));
                let (txx, txy, txz, fe) = visc_flux_x(
                    &grads[a], &grads[b], u[a], vv[a], ww[a], u[b], vv[b], ww[b], mu, kappa,
                );
                let (dxi, dmv, dms, de) = flux_deltas(txx, txy, txz, fe, dt / g.dx);
                apply_x(state, a, b, dxi, dmv, dms, de);
            }
            // wall faces at i=0 / i=nx when solid: one-sided wall shear.
            if walls.west {
                // low-side wall: the cell update is -f(wall) in the
                // divergence form, so the wall drains every component.
                let (a, b) = (idx(0, j, k), idx(1, j, k));
                let txx = 2.0 * mu * wall_grad(u[a], u[b], g.dx);
                let txy = mu * wall_grad(vv[a], vv[b], g.dx);
                let txz = mu * wall_grad(ww[a], ww[b], g.dx);
                let (dxi, dmv, dms, de) = flux_deltas(txx, txy, txz, 0.0, dt / g.dx);
                apply_x_wall_low(state, a, dxi, dmv, dms, de);
            }
            if walls.east {
                // high-side wall: the +x-oriented tangential stresses are
                // negative (the wall drags against the flow), and the cell
                // update is +f(wall), draining the cell.
                let (a, b) = (idx(nx - 1, j, k), idx(nx - 2, j, k));
                let txx = -2.0 * mu * wall_grad(u[a], u[b], g.dx);
                let txy = -mu * wall_grad(vv[a], vv[b], g.dx);
                let txz = -mu * wall_grad(ww[a], ww[b], g.dx);
                let (dxi, dmv, dms, de) = flux_deltas(txx, txy, txz, 0.0, dt / g.dx);
                apply_x_wall_high(state, a, dxi, dmv, dms, de);
            }
        }
    }
    // interior y-face fluxes: faces 1..ny-1 for every (i, k).
    for k in 0..nz {
        for j in 1..ny {
            for i in 0..nx {
                let (a, b) = (idx(i, j - 1, k), idx(i, j, k));
                let (tyx, tyy, tyz, fe) = visc_flux_y(
                    &grads[a], &grads[b], u[a], vv[a], ww[a], u[b], vv[b], ww[b], mu, kappa,
                );
                let (dxi, dmv, dms, de) = flux_deltas(tyx, tyy, tyz, fe, dt / g.dy);
                apply_y(state, a, b, dxi, dmv, dms, de);
            }
        }
    }
    // interior z-face fluxes: faces 1..nz-1 for every (i, j).
    for k in 1..nz {
        for j in 0..ny {
            for i in 0..nx {
                let (a, b) = ((k - 1) * ny * nx + j * nx + i, k * ny * nx + j * nx + i);
                let (tzx, tzy, tzz, fe) = visc_flux_z(
                    &grads[a], &grads[b], u[a], vv[a], ww[a], u[b], vv[b], ww[b], mu, kappa,
                );
                let (dxi, dmv, dms, de) = flux_deltas(tzx, tzy, tzz, fe, dt / g.dz);
                apply_z(state, a, b, dxi, dmv, dms, de);
            }
        }
    }
    // solid wall faces in y and z: one-sided shear, adiabatic (no heat).
    if walls.south {
        // low-side wall: cell -= f, tangential stresses positive.
        for k in 0..nz {
            for i in 0..nx {
                let (a, b) = (idx(i, 0, k), idx(i, 1, k));
                let (tyx, tyy, tyz) = (
                    mu * wall_grad(u[a], u[b], g.dy),
                    2.0 * mu * wall_grad(vv[a], vv[b], g.dy),
                    mu * wall_grad(ww[a], ww[b], g.dy),
                );
                let (dxi, dmv, dms, de) = flux_deltas(tyx, tyy, tyz, 0.0, dt / g.dy);
                apply_y_wall_low(state, a, dxi, dmv, dms, de);
            }
        }
    }
    if walls.north {
        // high-side wall: +y stresses negative, cell += f.
        for k in 0..nz {
            for i in 0..nx {
                let (a, b) = (idx(i, ny - 1, k), idx(i, ny - 2, k));
                let (tyx, tyy, tyz) = (
                    -mu * wall_grad(u[a], u[b], g.dy),
                    -2.0 * mu * wall_grad(vv[a], vv[b], g.dy),
                    -mu * wall_grad(ww[a], ww[b], g.dy),
                );
                let (dxi, dmv, dms, de) = flux_deltas(tyx, tyy, tyz, 0.0, dt / g.dy);
                apply_y_wall_high(state, a, dxi, dmv, dms, de);
            }
        }
    }
    if walls.bottom {
        // low-side wall: cell -= f, tangential stresses positive.
        for j in 0..ny {
            for i in 0..nx {
                let (a, b) = (idx(i, j, 0), idx(i, j, 1));
                let (tzx, tzy, tzz) = (
                    mu * wall_grad(u[a], u[b], g.dz),
                    mu * wall_grad(vv[a], vv[b], g.dz),
                    2.0 * mu * wall_grad(ww[a], ww[b], g.dz),
                );
                let (dxi, dmv, dms, de) = flux_deltas(tzx, tzy, tzz, 0.0, dt / g.dz);
                apply_z_wall_low(state, a, dxi, dmv, dms, de);
            }
        }
    }
    if walls.top {
        // high-side wall: +z stresses negative, cell += f.
        for j in 0..ny {
            for i in 0..nx {
                let (a, b) = (idx(i, j, nz - 1), idx(i, j, nz - 2));
                let (tzx, tzy, tzz) = (
                    -mu * wall_grad(u[a], u[b], g.dz),
                    -mu * wall_grad(vv[a], vv[b], g.dz),
                    -2.0 * mu * wall_grad(ww[a], ww[b], g.dz),
                );
                let (dxi, dmv, dms, de) = flux_deltas(tzx, tzy, tzz, 0.0, dt / g.dz);
                apply_z_wall_high(state, a, dxi, dmv, dms, de);
            }
        }
    }
    Ok(())
}

/// the y-face viscous flux: the tau_y stress row dotted with the face
/// velocity plus heat conduction.
#[allow(clippy::too_many_arguments)]
fn visc_flux_y(
    ga: &Grads3d,
    gb: &Grads3d,
    ua: f64,
    va: f64,
    wa: f64,
    ub: f64,
    vb: f64,
    wb: f64,
    mu: f64,
    kappa: f64,
) -> (f64, f64, f64, f64) {
    let (uy, vy, wy, ty) = (
        0.5 * (ga.uy + gb.uy),
        0.5 * (ga.vy + gb.vy),
        0.5 * (ga.wy + gb.wy),
        0.5 * (ga.ty + gb.ty),
    );
    let (ux, wz) = (0.5 * (ga.ux + gb.ux), 0.5 * (ga.wz + gb.wz));
    let (vx, vz) = (0.5 * (ga.vx + gb.vx), 0.5 * (ga.vz + gb.vz));
    let (uf, vf, wf) = (0.5 * (ua + ub), 0.5 * (va + vb), 0.5 * (wa + wb));
    let divu = ux + vy + wz;
    let tyy = mu * (2.0 * vy - 2.0 / 3.0 * divu);
    let tyx = mu * (uy + vx);
    let tyz = mu * (wy + vz);
    let q = kappa * ty;
    let fe = uf * tyx + vf * tyy + wf * tyz + q;
    (tyx, tyy, tyz, fe)
}

/// the z-face viscous flux: the tau_z stress row dotted with the face
/// velocity plus heat conduction.
#[allow(clippy::too_many_arguments)]
fn visc_flux_z(
    ga: &Grads3d,
    gb: &Grads3d,
    ua: f64,
    va: f64,
    wa: f64,
    ub: f64,
    vb: f64,
    wb: f64,
    mu: f64,
    kappa: f64,
) -> (f64, f64, f64, f64) {
    let (uz, vz, wz, tz) = (
        0.5 * (ga.uz + gb.uz),
        0.5 * (ga.vz + gb.vz),
        0.5 * (ga.wz + gb.wz),
        0.5 * (ga.tz + gb.tz),
    );
    let (ux, vy) = (0.5 * (ga.ux + gb.ux), 0.5 * (ga.vy + gb.vy));
    let (wx, wy) = (0.5 * (ga.wx + gb.wx), 0.5 * (ga.wy + gb.wy));
    let (uf, vf, wf) = (0.5 * (ua + ub), 0.5 * (va + vb), 0.5 * (wa + wb));
    let divu = ux + vy + wz;
    let tzz = mu * (2.0 * wz - 2.0 / 3.0 * divu);
    let tzx = mu * (uz + wx);
    let tzy = mu * (vz + wy);
    let q = kappa * tz;
    let fe = uf * tzx + vf * tzy + wf * tzz + q;
    (tzx, tzy, tzz, fe)
}

/// the per-axis flux difference contribution for one face, scaled by dt/d.
fn flux_deltas(f1: f64, f2: f64, f3: f64, fe: f64, s: f64) -> (f64, f64, f64, f64) {
    (s * f1, s * f2, s * f3, s * fe)
}

fn apply_y(st: &mut ConservedState3d, a: usize, b: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] += dxi;
    st.my[a] += dmv;
    st.mz[a] += dms;
    st.e[a] += de;
    st.mx[b] -= dxi;
    st.my[b] -= dmv;
    st.mz[b] -= dms;
    st.e[b] -= de;
}

fn apply_z(st: &mut ConservedState3d, a: usize, b: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] += dxi;
    st.my[a] += dmv;
    st.mz[a] += dms;
    st.e[a] += de;
    st.mx[b] -= dxi;
    st.my[b] -= dmv;
    st.mz[b] -= dms;
    st.e[b] -= de;
}

fn apply_y_wall_low(st: &mut ConservedState3d, a: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] -= dxi;
    st.my[a] -= dmv;
    st.mz[a] -= dms;
    st.e[a] -= de;
}

fn apply_y_wall_high(st: &mut ConservedState3d, a: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] += dxi;
    st.my[a] += dmv;
    st.mz[a] += dms;
    st.e[a] += de;
}

fn apply_z_wall_low(st: &mut ConservedState3d, a: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] -= dxi;
    st.my[a] -= dmv;
    st.mz[a] -= dms;
    st.e[a] -= de;
}

fn apply_z_wall_high(st: &mut ConservedState3d, a: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] += dxi;
    st.my[a] += dmv;
    st.mz[a] += dms;
    st.e[a] += de;
}

fn apply_x(st: &mut ConservedState3d, a: usize, b: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] += dxi;
    st.my[a] += dmv;
    st.mz[a] += dms;
    st.e[a] += de;
    st.mx[b] -= dxi;
    st.my[b] -= dmv;
    st.mz[b] -= dms;
    st.e[b] -= de;
}

fn apply_x_wall_low(st: &mut ConservedState3d, a: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] -= dxi;
    st.my[a] -= dmv;
    st.mz[a] -= dms;
    st.e[a] -= de;
}

fn apply_x_wall_high(st: &mut ConservedState3d, a: usize, dxi: f64, dmv: f64, dms: f64, de: f64) {
    st.mx[a] += dxi;
    st.my[a] += dmv;
    st.mz[a] += dms;
    st.e[a] += de;
}

/// second-order heun step for the 3d navier-stokes operator: inviscid
/// predictor stages with the viscous add at each stage, averaged at the end
/// (the 3d mirror of advance2d_visc_rk2).
pub fn advance3d_visc_rk2(
    state: &mut ConservedState3d,
    g: &Grid3d,
    gamma: f64,
    cfl: f64,
    muscl: bool,
    v: ViscParams3d,
    walls: &Walls3d,
) -> Result<(f64, f64), crate::Error> {
    let mu = v.mu;
    let pr = v.pr;
    let kappa = mu * gamma / ((gamma - 1.0) * pr);
    let rho_min = state
        .rho
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .max(1e-6);
    // the explicit diffusion bound caps the step.
    let d = g.dx.min(g.dy).min(g.dz);
    let dt_visc = 0.25 * d * d * rho_min / (mu * kappa.max(1.0));
    let mut s1 = state.clone();
    let (dt1, _) = {
        // the inviscid step capped by the viscous bound: emulate the 2d
        // capped path by scaling cfl so the resulting dt respects dt_visc.
        let (u, vv, ww, et) = cons_to_prim3d(&s1.rho, &s1.mx, &s1.my, &s1.mz, &s1.e)?;
        let p = eos_pressure3d(gamma, &s1.rho, &et, &u, &vv, &ww)?;
        let smax = smax3d_of(&u, &vv, &ww, &p, &s1.rho, gamma);
        let dt_acoustic = cfl * d / smax;
        let dt = dt_acoustic.min(dt_visc);
        march_with_dt(&mut s1, g, gamma, dt, muscl)?;
        add_viscous3d(&mut s1, g, gamma, v, walls, dt)?;
        (dt, 0.0)
    };
    let (dt2, _) = {
        let (u, vv, ww, et) = cons_to_prim3d(&s1.rho, &s1.mx, &s1.my, &s1.mz, &s1.e)?;
        let p = eos_pressure3d(gamma, &s1.rho, &et, &u, &vv, &ww)?;
        let smax = smax3d_of(&u, &vv, &ww, &p, &s1.rho, gamma);
        let dt_acoustic = cfl * d / smax;
        let dt = dt_acoustic.min(dt_visc);
        march_with_dt(&mut s1, g, gamma, dt, muscl)?;
        add_viscous3d(&mut s1, g, gamma, v, walls, dt)?;
        (dt, 0.0)
    };
    let _ = dt2;
    let n = g.nx * g.ny * g.nz;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.mz[k] = 0.5 * state.mz[k] + 0.5 * s1.mz[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt1, 0.0))
}

/// the largest fast-characteristic speed over the cells.
fn smax3d_of(u: &[f64], v: &[f64], w: &[f64], p: &[f64], rho: &[f64], gamma: f64) -> f64 {
    let mut s = 0.0f64;
    for k in 0..rho.len() {
        let a = (gamma * p[k].max(1e-12) / rho[k]).sqrt();
        s = s.max(u[k].abs().max(v[k].abs().max(w[k].abs())) + a);
    }
    s
}

/// one inviscid euler step with an externally fixed dt (the capped march
/// the rk2 viscous driver needs; solver3d's own advance3d computes its dt
/// internally, so this scales the cfl to land on the requested dt).
fn march_with_dt(
    state: &mut ConservedState3d,
    g: &Grid3d,
    gamma: f64,
    dt: f64,
    muscl: bool,
) -> Result<(), crate::Error> {
    let d = g.dx.min(g.dy).min(g.dz);
    let (u, v, w, et) = cons_to_prim3d(&state.rho, &state.mx, &state.my, &state.mz, &state.e)?;
    let p = eos_pressure3d(gamma, &state.rho, &et, &u, &v, &w)?;
    let smax = smax3d_of(&u, &v, &w, &p, &state.rho, gamma);
    let cfl = dt * smax / d;
    advance3d(state, g, gamma, cfl, muscl)?;
    Ok(())
}

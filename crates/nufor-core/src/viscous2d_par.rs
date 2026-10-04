//! the threaded sa march: the same heun-coupled mean-flow + transport step
//! as viscous2d::advance2d_sa_rk2, with the three heavy kernels sharded
//! across scoped threads and their deltas applied serially so the result is
//! bit-identical to the serial march at every thread count.

use crate::grid2d::Grid2d;
use crate::sa::eddy_viscosity;
use crate::solver2d::{advance2d_capped_dts, TimeControl};
use crate::turb2d::TurbState;
use crate::viscous2d::{add_viscous_cells, sa_dt_cap, TurbCtx};
use crate::{Boundaries2d, ConservedState2d, Error, ViscParams};

/// shard `n` independent work items across the pool; each thread writes its
/// own slab, collected in order so results never interleave.
fn shard<T: Send>(n: usize, nthreads: usize, work: impl Fn(usize) -> T + Sync) -> Vec<T> {
    let threads = nthreads.min(n).max(1);
    let chunk = n.div_ceil(threads);
    let mut out: Vec<Option<T>> = (0..n).map(|_| None).collect();
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for start in (0..n).step_by(chunk) {
            let end = (start + chunk).min(n);
            let work = &work;
            handles.push(s.spawn(move || (start, (start..end).map(work).collect::<Vec<_>>())));
        }
        for h in handles {
            let (start, v) = h.join().unwrap();
            for (k, item) in v.into_iter().enumerate() {
                out[start + k] = Some(item);
            }
        }
    });
    out.into_iter().map(|x| x.unwrap()).collect()
}

/// the step configuration shared by the threaded marches.
#[derive(Debug, Clone, Copy)]
pub struct StepConfig {
    pub gamma: f64,
    pub cfl: f64,
    pub muscl: bool,
    pub nthreads: usize,
}

/// the threaded sa heun step. `nthreads` <= 1 falls back to the serial march
/// so small grids never pay thread overhead.
pub fn advance2d_sa_rk2_par(
    state: &mut ConservedState2d,
    turb: &mut TurbState,
    g: &Grid2d,
    bc: &Boundaries2d,
    cfg: StepConfig,
) -> Result<f64, Error> {
    let (gamma, cfl, muscl, nthreads) = (cfg.gamma, cfg.cfl, cfg.muscl, cfg.nthreads);
    if nthreads <= 1 {
        return crate::viscous2d::advance2d_sa_rk2(state, turb, g, gamma, cfl, muscl, bc);
    }
    let n = g.nx * g.ny;
    let sa = turb.params;
    let d = &turb.d;
    if turb.nu_tilde.len() != n || d.len() != n {
        return Err(Error::InvalidArgs);
    }
    let nu = sa.mu;
    let nt_max = turb
        .nu_tilde
        .iter()
        .cloned()
        .fold(sa.nu_tilde_inf, f64::max);
    let nu_max = nu / state.rho.iter().cloned().fold(f64::INFINITY, f64::min) + nt_max;
    let dt_cap = sa_dt_cap(g, nu_max);
    let mu_t_of = |t: &TurbState, rho: &[f64]| -> Vec<f64> {
        t.nu_tilde
            .iter()
            .zip(rho)
            .map(|(&nt, &r)| eddy_viscosity(r, nt, nu / r.max(1e-12)))
            .collect()
    };
    // stage 1.
    let mut s1 = state.clone();
    let mut t1 = turb.clone();
    let (dt1, _) = crate::solver2d::advance2d_capped_par(
        &mut s1,
        g,
        gamma,
        cfl,
        muscl,
        bc,
        TimeControl::Global(dt_cap),
        nthreads,
    )?;
    let mu_t1 = mu_t_of(&t1, &s1.rho);
    let ctx1 = TurbCtx {
        mu_t: &mu_t1,
        pr_t: sa.pr_t,
        bc,
    };
    add_viscous_cells_par(
        &mut s1,
        g,
        gamma,
        nu,
        sa.pr,
        Some(ctx1),
        TimeControl::Global(dt1),
        nthreads,
    )?;
    advance_turb_par(&mut t1, &s1, g, bc, dt1, None, nthreads)?;
    // stage 2 from the stage-1 state.
    let (dt2, _) = crate::solver2d::advance2d_capped_par(
        &mut s1,
        g,
        gamma,
        cfl,
        muscl,
        bc,
        TimeControl::Global(dt_cap),
        nthreads,
    )?;
    let mu_t1b = mu_t_of(&t1, &s1.rho);
    let ctx2 = TurbCtx {
        mu_t: &mu_t1b,
        pr_t: sa.pr_t,
        bc,
    };
    add_viscous_cells_par(
        &mut s1,
        g,
        gamma,
        nu,
        sa.pr,
        Some(ctx2),
        TimeControl::Global(dt2),
        nthreads,
    )?;
    advance_turb_par(&mut t1, &s1, g, bc, dt2, None, nthreads)?;
    // heun average.
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
        turb.nu_tilde[k] = 0.5 * turb.nu_tilde[k] + 0.5 * t1.nu_tilde[k];
    }
    Ok(dt1)
}

/// the threaded laminar viscous heun step (the same shape as
/// advance2d_visc_rk2, sharded where it pays).
pub fn advance2d_visc_rk2_par(
    state: &mut ConservedState2d,
    g: &Grid2d,
    bc: &Boundaries2d,
    v: ViscParams,
    cfg: StepConfig,
) -> Result<(f64, f64), Error> {
    let (gamma, cfl, muscl, nthreads) = (cfg.gamma, cfg.cfl, cfg.muscl, cfg.nthreads);
    if nthreads <= 1 {
        return crate::viscous2d::advance2d_visc_rk2(state, g, gamma, cfl, muscl, bc, v);
    }
    let mu = v.mu;
    let pr = v.pr;
    let kappa = gamma / ((gamma - 1.0) * pr);
    let rho_min = state
        .rho
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .max(1e-6);
    let dt_visc = 0.25 * g.dx.min(g.dy).powi(2) * rho_min / (mu * kappa.max(1.0));
    let mut s1 = state.clone();
    let (dt1, _) = advance2d_capped_dts(
        &mut s1,
        g,
        gamma,
        cfl,
        muscl,
        bc,
        TimeControl::Global(dt_visc),
    )?;
    crate::viscous2d::add_viscous_bc(&mut s1, g, gamma, mu, pr, bc, dt1)?;
    let (dt2, _) = advance2d_capped_dts(
        &mut s1,
        g,
        gamma,
        cfl,
        muscl,
        bc,
        TimeControl::Global(dt_visc),
    )?;
    crate::viscous2d::add_viscous_bc(&mut s1, g, gamma, mu, pr, bc, dt2)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt1, 0.0))
}

/// the threaded sa transport: per-cell right-hand sides computed in
/// parallel, applied in cell order so the update is bit-identical.
pub fn advance_turb_par(
    turb: &mut TurbState,
    state: &ConservedState2d,
    g: &Grid2d,
    bc: &Boundaries2d,
    dt: f64,
    dts: Option<&[f64]>,
    nthreads: usize,
) -> Result<(), Error> {
    use crate::sa::{source, SIGMA};
    use crate::turb2d::vorticity;
    use crate::C_B2;
    let mu_lam = turb.params.mu;
    let d = &turb.d;
    let nx = g.nx;
    let ny = g.ny;
    let n = nx * ny;
    if turb.nu_tilde.len() != n || d.len() != n {
        return Err(Error::InvalidArgs);
    }
    let (u, v, _et) = crate::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let unif_x = (0..nx).all(|i| (g.dxs[i] - g.dxs[0]).abs() <= 1e-9 * g.dxs[0].abs().max(1e-12));
    let unif_y = (0..ny).all(|j| (g.dys[j] - g.dys[0]).abs() <= 1e-9 * g.dys[0].abs().max(1e-12));
    let s = vorticity(&u, &v, g, unif_x, unif_y);
    let pn = crate::turb2d::pad_scalar(&turb.nu_tilde, nx, ny, bc);
    let w = nx + 2;
    let inv_xs: Vec<f64> = g.dxs.iter().map(|wv| 1.0 / wv).collect();
    let inv_ys: Vec<f64> = g.dys.iter().map(|wv| 1.0 / wv).collect();
    let (inv_x, inv_y) = if unif_x && unif_y {
        (Some(1.0 / g.dx), Some(1.0 / g.dy))
    } else {
        (None, None)
    };
    // per-cell rhs values, computed in parallel.
    let rhs: Vec<f64> = shard(n, nthreads, |c| {
        let i = c % nx;
        let j = c / nx;
        let rho = state.rho[c].max(1e-12);
        let nu_c = mu_lam / rho;
        let mut f_adv = [0.0f64; 4];
        let mut f_dif = [0.0f64; 4];
        for (k, f) in [i, i + 1].iter().enumerate() {
            let f = *f;
            let (la, lb) = (f, f + 1);
            let (cell_l, cell_r) = (f as isize - 1, f as isize);
            let u_f = if f == 0 {
                u[j * nx]
            } else if f == nx {
                u[j * nx + nx - 1]
            } else {
                0.5 * (u[j * nx + cell_l as usize] + u[j * nx + cell_r as usize])
            };
            let rho_f = if f == 0 {
                state.rho[j * nx]
            } else if f == nx {
                state.rho[j * nx + nx - 1]
            } else {
                0.5 * (state.rho[j * nx + cell_l as usize] + state.rho[j * nx + cell_r as usize])
            };
            let left = pn[(j + 1) * w + la];
            let right = pn[(j + 1) * w + lb];
            let nt_up = if u_f >= 0.0 { left } else { right };
            f_adv[k] = rho_f * u_f * nt_up;
            let nt_f = 0.5 * (left + right);
            let nu_f = mu_lam / rho_f.max(1e-12);
            let inv = if let Some(ix) = inv_x {
                ix
            } else {
                let dc = if f == 0 {
                    g.centers_x[c] - g.faces_x[0]
                } else if f == nx {
                    g.faces_x[nx] - g.centers_x[j * nx + nx - 1]
                } else {
                    g.centers_x[j * nx + f] - g.centers_x[j * nx + f - 1]
                };
                1.0 / dc
            };
            f_dif[k] = (nu_f + nt_f) / SIGMA * (right - left) * inv;
        }
        for (k, f) in [j, j + 1].iter().enumerate() {
            let f = *f;
            let (ra, rb) = (f, f + 1);
            let (cell_b, cell_t) = (f as isize - 1, f as isize);
            let v_f = if f == 0 {
                v[i]
            } else if f == ny {
                v[(ny - 1) * nx + i]
            } else {
                0.5 * (v[cell_b as usize * nx + i] + v[cell_t as usize * nx + i])
            };
            let rho_f = if f == 0 {
                state.rho[i]
            } else if f == ny {
                state.rho[(ny - 1) * nx + i]
            } else {
                0.5 * (state.rho[cell_b as usize * nx + i] + state.rho[cell_t as usize * nx + i])
            };
            let below = pn[ra * w + i + 1];
            let above = pn[rb * w + i + 1];
            let nt_up = if v_f >= 0.0 { below } else { above };
            f_adv[k + 2] = rho_f * v_f * nt_up;
            let nt_f = 0.5 * (below + above);
            let nu_f = mu_lam / rho_f.max(1e-12);
            let inv = if let Some(iy) = inv_y {
                iy
            } else {
                let dc = if f == 0 {
                    g.centers_y[c] - g.faces_y[0]
                } else if f == ny {
                    g.faces_y[ny] - g.centers_y[(ny - 1) * nx + i]
                } else {
                    g.centers_y[f * nx + i] - g.centers_y[(f - 1) * nx + i]
                };
                1.0 / dc
            };
            f_dif[k + 2] = (nu_f + nt_f) / SIGMA * (above - below) * inv;
        }
        let inv_x_i = if unif_x {
            inv_x.unwrap_or(1.0 / g.dx)
        } else {
            inv_xs[i]
        };
        let inv_y_j = if unif_y {
            inv_y.unwrap_or(1.0 / g.dy)
        } else {
            inv_ys[j]
        };
        let adv_net = (f_adv[1] - f_adv[0]) * inv_x_i + (f_adv[3] - f_adv[2]) * inv_y_j;
        let dif_net = (f_dif[1] - f_dif[0]) * inv_x_i + (f_dif[3] - f_dif[2]) * inv_y_j;
        let (ip, jp) = (i + 1, j + 1);
        let (gx, gy) = if unif_x && unif_y {
            (
                (pn[jp * w + i + 2] - pn[jp * w + i]) / (2.0 * g.dx),
                (pn[(j + 2) * w + ip] - pn[j * w + ip]) / (2.0 * g.dy),
            )
        } else {
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
                    pn[jp * w + i + 2],
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
        -adv_net / rho + dif_net + cb2_term + src
    });
    for c in 0..n {
        let dt_c = dts.map(|a| a[c]).unwrap_or(dt);
        turb.nu_tilde[c] = (turb.nu_tilde[c] + dt_c * rhs[c]).max(0.0);
    }
    Ok(())
}

/// the threaded viscous add: x-face rows and y-face columns sharded across
/// threads into the same flux arrays the serial operator builds, then the
/// per-cell divergence applied in cell order (bit-identical).
#[allow(clippy::too_many_arguments)]
pub fn add_viscous_cells_par(
    state: &mut ConservedState2d,
    g: &Grid2d,
    gamma: f64,
    mu_lam: f64,
    pr: f64,
    turb: Option<TurbCtx>,
    tc: TimeControl,
    nthreads: usize,
) -> Result<(), Error> {
    use crate::viscous2d::{
        cell_grads, cell_grads_metric, face_flux_split, face_grad_interp, pad2d, pad2d_vel,
        wall_cell_grad, wall_cell_grad_st, wall_shear_deriv, wall_shear_flux, MeshMetrics,
        ViscCoeffs, LAMINAR_BC,
    };
    // measured on the clustered plate: thread-spawn overhead swamps the
    // parallel gain below ~20k cells, so small meshes stay serial.
    if nthreads <= 1 || g.nx * g.ny < 20_000 {
        return add_viscous_cells(state, g, gamma, mu_lam, pr, turb, tc);
    }
    let (dt, dts) = tc.parts();
    let bc = turb.map(|t| t.bc).unwrap_or(&LAMINAR_BC);
    let mu_t: &[f64] = turb.map(|t| t.mu_t).unwrap_or(&[]);
    let pr_t = turb.map(|t| t.pr_t).unwrap_or(pr);
    let nx = g.nx;
    let ny = g.ny;
    let idx = |i: usize, j: usize| j * nx + i;
    let is_turb = !mu_t.is_empty();
    if is_turb && mu_t.len() != nx * ny {
        return Err(Error::InvalidArgs);
    }
    if dt <= 0.0 {
        return Ok(());
    }
    let (u, v, et) = crate::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let p = crate::eos_pressure2d(gamma, &state.rho, &et, &u, &v)?;
    let t: Vec<f64> = p.iter().zip(&state.rho).map(|(pp, r)| pp / r).collect();
    let centers_y: Vec<f64> = (0..ny).map(|j| g.centers_y[j * nx]).collect();
    let (pu, pv, pt) = (
        pad2d_vel(&u, nx, ny, bc, &centers_y, 0),
        pad2d_vel(&v, nx, ny, bc, &centers_y, 1),
        pad2d(&t, nx, ny),
    );
    let unif_x = (0..nx).all(|i| (g.dxs[i] - g.dxs[0]).abs() <= 1e-9 * g.dxs[0].abs().max(1e-12));
    let unif_y = (0..ny).all(|j| (g.dys[j] - g.dys[0]).abs() <= 1e-9 * g.dys[0].abs().max(1e-12));
    let m = if unif_x && unif_y {
        None
    } else {
        Some(MeshMetrics::new(g))
    };
    // gradients: rows are independent, shard them.
    let mut grads = match &m {
        None => cell_grads(&pu, &pv, &pt, nx, ny, 2.0 * g.dx, 2.0 * g.dy),
        Some(mm) => cell_grads_metric(&pu, &pv, &pt, nx, ny, mm),
    };
    let south = matches!(bc.south, crate::solver2d::Bc2d::NoSlipWall);
    let north = matches!(bc.north, crate::solver2d::Bc2d::NoSlipWall);
    if south {
        for i in 0..nx {
            let (a, b) = (idx(i, 0), idx(i, 1));
            match &m {
                None => wall_cell_grad(&mut grads[a], u[a], u[b], v[a], v[b], g.dy, 1.0),
                Some(mm) => wall_cell_grad_st(
                    &mut grads[a],
                    (u[a], v[a]),
                    (u[b], v[b]),
                    g.faces_y[0],
                    mm.y_c[0],
                    mm.y_c[1],
                ),
            }
        }
    }
    if north {
        for i in 0..nx {
            let (a, b) = (idx(i, ny - 1), idx(i, ny - 2));
            match &m {
                None => wall_cell_grad(&mut grads[a], u[a], u[b], v[a], v[b], g.dy, -1.0),
                Some(mm) => wall_cell_grad_st(
                    &mut grads[a],
                    (u[a], v[a]),
                    (u[b], v[b]),
                    g.faces_y[ny],
                    mm.y_c[ny - 1],
                    mm.y_c[ny - 2],
                ),
            }
        }
    }
    let kappa = gamma / ((gamma - 1.0) * pr);
    let kappa_t = if is_turb {
        gamma / ((gamma - 1.0) * pr_t)
    } else {
        0.0
    };
    let coeffs = ViscCoeffs {
        kappa,
        kappa_t,
        gamma,
    };
    let (dtdx, dtdy) = (dt / g.dx, dt / g.dy);
    // x-face fluxes: one row per j, interior faces only.
    let fxm_fym_fe: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = shard(ny, nthreads, |j| {
        let mut fxm = vec![0.0; nx + 1];
        let mut fym = vec![0.0; nx + 1];
        let mut fe = vec![0.0; nx + 1];
        for f in 1..nx {
            let (a, b) = (idx(f - 1, j), idx(f, j));
            let gf = face_grad_interp(
                &grads[a],
                &grads[b],
                u[a],
                u[b],
                v[a],
                v[b],
                (g.centers_x[a], g.centers_x[b], g.faces_x[f]),
            );
            let (txx, txy, _, qx, _) = face_flux_split(&gf, mu_lam, mu_t, a, b, &coeffs);
            fxm[f] = txx;
            fym[f] = txy;
            fe[f] = gf.u * txx + gf.v * txy + qx;
        }
        (fxm, fym, fe)
    });
    let mut fxm = vec![0.0; (nx + 1) * ny];
    let mut fym = vec![0.0; (nx + 1) * ny];
    let mut fe = vec![0.0; (nx + 1) * ny];
    for (j, row) in fxm_fym_fe.iter().enumerate() {
        let (a, b, c) = (&row.0, &row.1, &row.2);
        for f in 1..nx {
            let k = f * ny + j;
            fxm[k] = a[f];
            fym[k] = b[f];
            fe[k] = c[f];
        }
    }
    // y-face fluxes: one column per i, including the wall faces.
    let gy: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = shard(nx, nthreads, |i| {
        let mut gym_x = vec![0.0; ny + 1];
        let mut gym_y = vec![0.0; ny + 1];
        let mut gye = vec![0.0; ny + 1];
        let y_c0 = g.centers_y[0];
        let y_c1 = g.centers_y[g.nx];
        let y_c_n = g.centers_y[(g.ny - 1) * g.nx];
        let y_c_nm1 = g.centers_y[(g.ny - 2) * g.nx];
        for f in 1..ny {
            let (a, b) = (idx(i, f - 1), idx(i, f));
            let cy_a = g.centers_y[a];
            let cy_b = g.centers_y[b];
            let gf = face_grad_interp(
                &grads[a],
                &grads[b],
                u[a],
                u[b],
                v[a],
                v[b],
                (cy_a, cy_b, g.faces_y[f]),
            );
            let (_, txy, tyy, _, qy) = face_flux_split(&gf, mu_lam, mu_t, a, b, &coeffs);
            gym_x[f] = txy;
            gym_y[f] = tyy;
            gye[f] = gf.u * txy + gf.v * tyy + qy;
        }
        if matches!(bc.south, crate::solver2d::Bc2d::NoSlipWall) {
            let (a, b) = (idx(i, 0), idx(i, 1));
            let mu_f = if mu_t.is_empty() {
                mu_lam
            } else {
                mu_lam + mu_t[a]
            };
            let (txy, tyy) = if unif_y {
                wall_shear_flux(u[a], u[b], v[a], v[b], mu_f, g.dy)
            } else {
                (
                    mu_f * wall_shear_deriv(u[a], u[b], g.faces_y[0], y_c0, y_c1),
                    mu_f * wall_shear_deriv(v[a], v[b], g.faces_y[0], y_c0, y_c1),
                )
            };
            gym_x[0] = txy;
            gym_y[0] = tyy;
            gye[0] = u[a] * txy + v[a] * tyy;
        }
        if matches!(bc.north, crate::solver2d::Bc2d::NoSlipWall) {
            let (a, b) = (idx(i, ny - 1), idx(i, ny - 2));
            let mu_f = if mu_t.is_empty() {
                mu_lam
            } else {
                mu_lam + mu_t[a]
            };
            let (txy, tyy) = if unif_y {
                wall_shear_flux(u[a], u[b], v[a], v[b], mu_f, g.dy)
            } else {
                (
                    -mu_f * wall_shear_deriv(u[a], u[b], g.faces_y[ny], y_c_n, y_c_nm1),
                    -mu_f * wall_shear_deriv(v[a], v[b], g.faces_y[ny], y_c_n, y_c_nm1),
                )
            };
            gym_x[ny] = -txy;
            gym_y[ny] = -tyy;
            gye[ny] = -(u[a] * txy + v[a] * tyy);
        }
        (gym_x, gym_y, gye)
    });
    let mut gym_x = vec![0.0; nx * (ny + 1)];
    let mut gym_y = vec![0.0; nx * (ny + 1)];
    let mut gye = vec![0.0; nx * (ny + 1)];
    for (i, col) in gy.iter().enumerate() {
        let (a, b, c) = (&col.0, &col.1, &col.2);
        for f in 0..=ny {
            let k = i * (ny + 1) + f;
            gym_x[k] = a[f];
            gym_y[k] = b[f];
            gye[k] = c[f];
        }
    }
    // the divergence apply, serially in cell order.
    for (j, &dyw) in g.dys.iter().enumerate() {
        for (i, &dxw) in g.dxs.iter().enumerate() {
            let c = idx(i, j);
            let (dxi, dyj) = match dts {
                None => (
                    if unif_x { dtdx } else { dt / dxw },
                    if unif_y { dtdy } else { dt / dyw },
                ),
                Some(a) => (a[c] / dxw, a[c] / dyw),
            };
            let (a, b) = (i * ny + j, (i + 1) * ny + j);
            state.mx[c] += dxi * (fxm[b] - fxm[a]);
            state.my[c] += dxi * (fym[b] - fym[a]);
            state.e[c] += dxi * (fe[b] - fe[a]);
            let (dd, ee) = (i * (ny + 1) + j, i * (ny + 1) + j + 1);
            state.mx[c] += dyj * (gym_x[ee] - gym_x[dd]);
            state.my[c] += dyj * (gym_y[ee] - gym_y[dd]);
            state.e[c] += dyj * (gye[ee] - gye[dd]);
        }
    }
    Ok(())
}

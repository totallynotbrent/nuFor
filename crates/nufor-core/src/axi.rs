//! the axisymmetric euler march: the 2d planar flux machinery with the
//! geometric source terms of a body of revolution.
//!
//! in cylindrical coordinates (x axial, r radial) the axisymmetric euler
//! equations keep the same face fluxes as the planar 2d case per unit
//! area, but the finite volume is an annulus, so the radial flux
//! difference carries the face radii and the radial momentum feels a
//! pressure source:
//!
//!   d/dt (rho) = -[ d(f_x)/dx + (1/r) d(r f_r)/dr ]
//!   d/dt (m_x) = -[ d(f_x^mx)/dx + (1/r) d(r f_r^mx)/dr ]
//!   d/dt (m_r) = -[ d(f_x^my)/dx + (1/r) d(r f_r^my)/dr ] + p / r
//!   d/dt (e)   = -[ d(f_x^e)/dx + (1/r) d(r f_r^e)/dr ]
//!
//! discretized on cell (i,j) with center radius r_c and radial face
//! radii r_- = r_j, r_+ = r_{j+1}:
//!
//!   u_t = -( f_x difference + (r_+ F_r[j+1] - r_- F_r[j]) / (r_c dr) )
//!         + (0, 0, p / r_c, 0)
//!
//! two consequences worth writing down: the south face of the domain is
//! the symmetry axis, and a slip-wall ghost (v -> -v, others copied)
//! makes the axis flux pure pressure transport with zero mass flow,
//! exactly the axis condition; and at the axis face r_- = 0, so the
//! first radial row only sees its outer face.
//!
//! the march reuses the planar flux machinery unchanged (fluxes per
//! unit area are identical in both formulations) and rewrites only the
//! update. the thermodynamic closure arrives as a `ThermoModel`: the
//! perfect-gas branch reconstructs pressure at faces (the historical
//! behavior), the equilibrium-air branch reconstructs internal energy
//! and evaluates the closure per face state.

use crate::grid2d::Grid2d;
use crate::hllc2d::{hllc_flux, FacePrim};
use crate::state2d::{cons_to_prim2d, ConservedState2d};
use crate::thermo::{self, ThermoModel};
use crate::Error;

use crate::solver2d::{ghost_value, Bc2d, Boundaries2d};

/// the van leer limiter: the harmonic mean, zero when the slopes disagree.
fn van_leer(dm: f64, dp: f64) -> f64 {
    if dm * dp <= 0.0 {
        0.0
    } else {
        2.0 * dm * dp / (dm + dp)
    }
}

/// muscl face states (left, right; each length n+1) from a padded slice
/// of length n+2, identical to the planar march's version.
fn face_states(p: &[f64], n: usize, muscl: bool) -> (Vec<f64>, Vec<f64>) {
    let mut d = vec![0.0; n];
    if muscl {
        for i in 0..n {
            let dm = p[i + 1] - p[i];
            let dp = p[i + 2] - p[i + 1];
            d[i] = van_leer(dm, dp);
        }
    }
    let mut fl = vec![0.0; n + 1];
    let mut fr = vec![0.0; n + 1];
    fl[0] = p[0];
    fr[0] = p[1] - 0.5 * d[0];
    fl[n] = p[n] + 0.5 * d[n - 1];
    fr[n] = p[n + 1];
    for f in 1..n {
        fl[f] = p[f] + 0.5 * d[f - 1];
        fr[f] = p[f + 1] - 0.5 * d[f];
    }
    (fl, fr)
}

/// the model-aware laminar viscous step: euler sub-steps through
/// advance2d_model_capped and the diffusive half-steps through
/// add_viscous_cells_mu with a per-cell sutherland viscosity from the
/// closure temperature. wall_temperature > 0 sets the isothermal
/// no-slip wall state; 0 means adiabatic.
pub fn advance2d_model_visc_rk2(
    state: &mut ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    mu0: f64,
    t0_ref: f64,
    s_param: f64,
    pr: f64,
    wall_temperature: f64,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let n = g.nx * g.ny;
    // the diffusion bound with the max kinematic viscosity the shock
    // layer reaches; using the reference-state mu underestimates the
    // stiffness at hot cells, so probe the actual field.
    let (u0, v0, et0) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0; n];
    for k in 0..n {
        e_int[k] = et0[k] - 0.5 * (u0[k] * u0[k] + v0[k] * v0[k]);
    }
    // fused (p, t, mu) pass; on eqair this replaces the temperature
    // bisection with one table read per cell.
    let (p0, _t0, mu_cells) = match model {
        ThermoModel::EqAir => crate::eqair_cea::fused_state_mu(&state.rho, &e_int, |t| {
            crate::viscous2d::sutherland_mu(mu0, t.max(1.0), t0_ref, s_param)
        }),
        ThermoModel::Perfect { .. } => {
            let p0 = thermo::pressure(model, &state.rho, &e_int, &u0, &v0)?;
            let t0 = thermo::temperature(model, &state.rho, &p0)?;
            let mut mc = vec![mu0; n];
            for k in 0..n {
                mc[k] = crate::viscous2d::sutherland_mu(mu0, t0[k].max(1.0), t0_ref, s_param);
            }
            (p0, t0, mc)
        }
    };
    let rho_min = state
        .rho
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .max(1e-6);
    let mu_max = mu_cells.iter().cloned().fold(0.0f64, f64::max);
    let nu_max = mu_max / rho_min;
    let dt_visc = 0.25 * g.dx.min(g.dy).powi(2) / nu_max.max(1e-12);
    let mut s1 = state.clone();
    let (dt1, _) =
        advance2d_model_capped(&mut s1, g, model, cfl, muscl, bc, dt_visc, nthreads, wall)?;
    crate::viscous2d::add_viscous_cells_mu(
        &mut s1,
        g,
        match model {
            ThermoModel::Perfect { gamma } => gamma,
            ThermoModel::EqAir => 1.4,
        },
        mu0,
        &mu_cells,
        pr,
        Some(crate::viscous2d::TurbCtx {
            mu_t: &[],
            pr_t: pr,
            bc,
        }),
        crate::solver2d::TimeControl::Global(dt1),
        Some(model),
        wall_temperature,
    )?;
    let (dt2, _) =
        advance2d_model_capped(&mut s1, g, model, cfl, muscl, bc, dt_visc, nthreads, wall)?;
    crate::viscous2d::add_viscous_cells_mu(
        &mut s1,
        g,
        match model {
            ThermoModel::Perfect { gamma } => gamma,
            ThermoModel::EqAir => 1.4,
        },
        mu0,
        &mu_cells,
        pr,
        Some(crate::viscous2d::TurbCtx {
            mu_t: &[],
            pr_t: pr,
            bc,
        }),
        crate::solver2d::TimeControl::Global(dt2),
        Some(model),
        wall_temperature,
    )?;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt1, 0.0))
}

/// closure-aware coupled spalart-allmaras march: the viscous mean-flow
/// advance above plus the sa transport, each heun stage consuming the
/// other's updated state. the eddy viscosity is rebuilt from nu_tilde
/// every stage so the momentum diffusion and the transported field stay
/// consistent. this is the eqair-capable twin of `advance2d_sa_rk2`,
/// which is hardcoded to a perfect-gas gamma.
pub fn advance2d_sa_model_rk2(
    state: &mut ConservedState2d,
    turb: &mut crate::turb2d::TurbState,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &crate::solver2d::Boundaries2d,
    mu0: f64,
    t0_ref: f64,
    s_param: f64,
    wall_temperature: f64,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<f64, Error> {
    use crate::sa::eddy_viscosity;
    use crate::turb2d::advance_turb;
    let n = g.nx * g.ny;
    let sa = turb.params;
    if turb.nu_tilde.len() != n || turb.d.len() != n {
        return Err(Error::InvalidArgs);
    }
    let gamma_eff = match model {
        ThermoModel::Perfect { gamma } => gamma,
        ThermoModel::EqAir => 1.4,
    };
    // eddy viscosity per cell from a nu_tilde field against a rho field.
    let mu_t_of = |nt: &[f64], rho: &[f64]| -> Vec<f64> {
        nt.iter()
            .zip(rho)
            .map(|(&ntk, &r)| eddy_viscosity(r, ntk, sa.mu / r.max(1e-12)))
            .collect()
    };
    // the step cap spans the stiffest of momentum diffusion and the sa
    // transport diffusivity over the (possibly hot) shock layer.
    let rho_min = state
        .rho
        .iter()
        .cloned()
        .fold(f64::INFINITY, f64::min)
        .max(1e-6);
    let nt_max = turb
        .nu_tilde
        .iter()
        .cloned()
        .fold(sa.nu_tilde_inf, f64::max);
    let nu_sa = sa.mu / rho_min + nt_max;
    let dt_sa = crate::viscous2d::sa_dt_cap(g, nu_sa);
    // stage 1: euler + viscous + turb from (state, turb).
    let mut s1 = state.clone();
    let mut t1 = turb.clone();
    let (u1, v1, et1) = cons_to_prim2d(&s1.rho, &s1.mx, &s1.my, &s1.e)?;
    let mut e1 = vec![0.0; n];
    for k in 0..n {
        e1[k] = et1[k] - 0.5 * (u1[k] * u1[k] + v1[k] * v1[k]);
    }
    let p1 = thermo::pressure(model, &s1.rho, &e1, &u1, &v1)?;
    let tt1 = thermo::temperature(model, &s1.rho, &p1)?;
    let mut mu_cells1 = vec![mu0; n];
    for k in 0..n {
        mu_cells1[k] = crate::viscous2d::sutherland_mu(mu0, tt1[k].max(1.0), t0_ref, s_param);
    }
    let (dt1, _) =
        advance2d_model_capped(&mut s1, g, model, cfl, muscl, bc, dt_sa, nthreads, wall)?;
    let mu_t1 = mu_t_of(&t1.nu_tilde, &s1.rho);
    crate::viscous2d::add_viscous_cells_mu(
        &mut s1,
        g,
        gamma_eff,
        mu0,
        &mu_cells1,
        sa.pr,
        Some(crate::viscous2d::TurbCtx {
            mu_t: &mu_t1,
            pr_t: sa.pr_t,
            bc,
        }),
        crate::solver2d::TimeControl::Global(dt1),
        Some(model),
        wall_temperature,
    )?;
    advance_turb(&mut t1, &s1, g, bc, dt1, None)?;
    // stage 2 re-evaluates both operators at the stage-1 state (heun).
    let (u2, v2, et2) = cons_to_prim2d(&s1.rho, &s1.mx, &s1.my, &s1.e)?;
    let mut e2 = vec![0.0; n];
    for k in 0..n {
        e2[k] = et2[k] - 0.5 * (u2[k] * u2[k] + v2[k] * v2[k]);
    }
    let p2 = thermo::pressure(model, &s1.rho, &e2, &u2, &v2)?;
    let tt2 = thermo::temperature(model, &s1.rho, &p2)?;
    let mut mu_cells2 = vec![mu0; n];
    for k in 0..n {
        mu_cells2[k] = crate::viscous2d::sutherland_mu(mu0, tt2[k].max(1.0), t0_ref, s_param);
    }
    let (dt2, _) =
        advance2d_model_capped(&mut s1, g, model, cfl, muscl, bc, dt_sa, nthreads, wall)?;
    let mu_t2 = mu_t_of(&t1.nu_tilde, &s1.rho);
    crate::viscous2d::add_viscous_cells_mu(
        &mut s1,
        g,
        gamma_eff,
        mu0,
        &mu_cells2,
        sa.pr,
        Some(crate::viscous2d::TurbCtx {
            mu_t: &mu_t2,
            pr_t: sa.pr_t,
            bc,
        }),
        crate::solver2d::TimeControl::Global(dt2),
        Some(model),
        wall_temperature,
    )?;
    advance_turb(&mut t1, &s1, g, bc, dt2, None)?;
    // heun average of the original and the twice-advanced state.
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
        turb.nu_tilde[k] = 0.5 * turb.nu_tilde[k] + 0.5 * t1.nu_tilde[k];
    }
    Ok(dt1)
}

/// like `advance2d_model` but with an external step cap, the viscous
/// path's way to interleave euler and diffusive half-steps at the
/// stiffer of the two bounds. single-advance.
fn advance2d_model_capped(
    state: &mut ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &crate::solver2d::Boundaries2d,
    dt_cap: f64,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let n = nx * ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0; n];
    // the prim->closure pass is per-cell: shard it with disjoint windows.
    crate::solver2d::shard_slice(&mut e_int, nthreads, |s, win| {
        for (w, k) in win.iter_mut().zip(s..) {
            *w = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        }
    });
    // one fused (p, a) table lookup per cell instead of two separate
    // log10-pair sweeps; both outputs are written through disjoint windows.
    let (mut p, mut a) = (vec![0.0f64; n], vec![0.0f64; n]);
    {
        let e_int = &e_int[..];
        let threads = if nthreads <= 1 || n < 2 {
            1
        } else {
            nthreads.min(n)
        };
        let chunk = n.div_ceil(threads);
        let mut a_rest = a.as_mut_slice();
        let mut p_rest = p.as_mut_slice();
        let mut offs: Vec<(usize, &mut [f64], &mut [f64])> = Vec::new();
        let mut start = 0usize;
        while start < n {
            let end = (start + chunk).min(n);
            let (ph, pt) = p_rest.split_at_mut(end - start);
            let (ah, at) = a_rest.split_at_mut(end - start);
            offs.push((start, ph, ah));
            p_rest = pt;
            a_rest = at;
            start = end;
        }
        let rho = &state.rho[..];
        std::thread::scope(|sc| {
            for (s, pw, aw) in offs {
                let e_int = &e_int[..];
                sc.spawn(move || {
                    for ((w, wa), k) in pw.iter_mut().zip(aw.iter_mut()).zip(s..) {
                        let (pk, ak) = match model {
                            ThermoModel::EqAir => crate::eqair_cea::p_a_at(rho[k], e_int[k]),
                            ThermoModel::Perfect { gamma } => {
                                let pk = (gamma - 1.0) * rho[k] * e_int[k];
                                let ak = (gamma * pk / rho[k].max(1e-12)).sqrt();
                                (pk, ak)
                            }
                        };
                        *w = pk;
                        *wa = ak;
                    }
                });
            }
        });
    }
    let smax = crate::solver2d::shard_rows_max(n, nthreads, |k| u[k].abs().max(v[k].abs()) + a[k]);
    let dt = (cfl * g.dx.min(g.dy) / smax.max(1e-12)).min(dt_cap);

    let (fx, fy) = compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, nthreads, wall)?;

    // the cell updates are disjoint per cell: shard each conserved array
    // through split windows and max-reduce the residual per chunk so the
    // result matches the serial loop exactly.
    let mut d_r = vec![0.0f64; n];
    let mut d_m = vec![0.0f64; n];
    let mut d_n = vec![0.0f64; n];
    let mut d_e = vec![0.0f64; n];
    let mut resids = vec![0.0f64; n];
    {
        let fxr = &fx[..];
        let fyr = &fy[..];
        crate::solver2d::shard_slice(&mut d_r, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let dxi = dt / g.dxs[i];
                let dyi = dt / g.dys[j];
                *w = -(dxi * (fxr[j].mass[i + 1] - fxr[j].mass[i])
                    + dyi * (fyr[i].mass[j + 1] - fyr[i].mass[j]));
            }
        });
        crate::solver2d::shard_slice(&mut d_m, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let dxi = dt / g.dxs[i];
                let dyi = dt / g.dys[j];
                *w = -(dxi * (fxr[j].mx[i + 1] - fxr[j].mx[i])
                    + dyi * (fyr[i].mx[j + 1] - fyr[i].mx[j]));
            }
        });
        crate::solver2d::shard_slice(&mut d_n, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let dxi = dt / g.dxs[i];
                let dyi = dt / g.dys[j];
                *w = -(dxi * (fxr[j].my[i + 1] - fxr[j].my[i])
                    + dyi * (fyr[i].my[j + 1] - fyr[i].my[j]));
            }
        });
        crate::solver2d::shard_slice(&mut d_e, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let dxi = dt / g.dxs[i];
                let dyi = dt / g.dys[j];
                *w = -(dxi * (fxr[j].e[i + 1] - fxr[j].e[i])
                    + dyi * (fyr[i].e[j + 1] - fyr[i].e[j]));
            }
        });
        let (dr, dm, dn, de) = (&d_r[..], &d_m[..], &d_n[..], &d_e[..]);
        crate::solver2d::shard_slice(&mut resids, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                *w = dr[k]
                    .abs()
                    .max(dm[k].abs())
                    .max(dn[k].abs())
                    .max(de[k].abs());
            }
        });
        // apply serially: exact serial-equivalent, cache-friendly, and the
        // deltas are already computed in parallel above.
        for k in 0..n {
            state.rho[k] += dr[k];
            state.mx[k] += dm[k];
            state.my[k] += dn[k];
            state.e[k] += de[k];
        }
    }
    let resid = resids.iter().cloned().fold(0.0f64, f64::max);
    Ok((dt, resid))
}

/// one axisymmetric euler step; compose with two half-steps for heun the
/// same way `advance2d_rk2` composes `advance2d`.
///
/// the south boundary is the symmetry axis: pass `Bc2d::SlipWall` there
/// (the ghost's negated v is exactly the axis condition).
pub fn advance2d_axi(
    state: &mut ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let n = nx * ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
    }
    let p = thermo::pressure(model, &state.rho, &e_int, &u, &v)?;
    let a = thermo::sound_speed(model, &state.rho, &et, &u, &v)?;
    let mut smax = 0.0f64;
    for k in 0..n {
        smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
    }
    let dt = cfl * g.dx.min(g.dy) / smax.max(1e-12);

    let (fx, fy) = compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, nthreads, wall)?;

    // the annular update: r-weighted radial flux difference plus the
    // radial pressure source.
    let mut resid = 0.0f64;
    for j in 0..ny {
        let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
        let rc = 0.5 * (rm + rp);
        let dr = g.dys[j];
        for (i, fycol) in fy.iter().enumerate() {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let radial = |f: &Vec<f64>| (rp * f[j + 1] - rm * f[j]) / (rc * dr);
            let drho = -(dxi * (fx[j].mass[i + 1] - fx[j].mass[i]) + dt * radial(&fycol.mass));
            let dmx = -(dxi * (fx[j].mx[i + 1] - fx[j].mx[i]) + dt * radial(&fycol.mx));
            let dmy =
                -(dxi * (fx[j].my[i + 1] - fx[j].my[i]) + dt * radial(&fycol.my)) + dt * p[k] / rc;
            let de = -(dxi * (fx[j].e[i + 1] - fx[j].e[i]) + dt * radial(&fycol.e));
            state.rho[k] += drho;
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            resid = resid
                .max(drho.abs())
                .max(dmx.abs())
                .max(dmy.abs())
                .max(de.abs());
        }
    }
    Ok((dt, resid))
}

/// face flux arrays for one sweep direction.
pub(crate) struct Sweep {
    pub(crate) mass: Vec<f64>,
    pub(crate) mx: Vec<f64>,
    pub(crate) my: Vec<f64>,
    pub(crate) e: Vec<f64>,
}

impl Sweep {
    fn new(n: usize) -> Self {
        Sweep {
            mass: vec![0.0; n + 1],
            mx: vec![0.0; n + 1],
            my: vec![0.0; n + 1],
            e: vec![0.0; n + 1],
        }
    }
}

/// build the hllc face state pair for one face from reconstructed
/// fields under the chosen closure. `x` is the reconstructed closure
/// variable: pressure for the perfect gas, internal energy per mass
/// for equilibrium air. returns (left, right).
#[allow(clippy::too_many_arguments)]
fn face_pair(
    model: ThermoModel,
    rl: f64,
    ul: f64,
    vl: f64,
    xl: f64,
    rr: f64,
    ur: f64,
    vr: f64,
    xr: f64,
) -> (FacePrim, FacePrim) {
    match model {
        ThermoModel::Perfect { gamma } => (
            FacePrim::perfect(rl, ul, vl, xl, gamma),
            FacePrim::perfect(rr, ur, vr, xr, gamma),
        ),
        ThermoModel::EqAir => {
            // the legacy tgas1 fits extrapolate past their window instead of
            // clamping, which the riemann solver needs: wake states exceed
            // the CEA table's energy ceiling and a clamped face pressure
            // detonates the base flow. consistency with the CEA state
            // advance is handled by the closure-variable reconstruction
            // (e_int) plus the fused cell pass.
            // the legacy tgas1 fits extrapolate past their window instead of
            // clamping, which the riemann solver needs: wake states exceed
            // the CEA table's energy ceiling and a clamped face pressure
            // detonates the base flow (re-verified: the wall zeroing does
            // not cure it — the offending faces are fluid-fluid in the
            // wake). consistency with the CEA state advance is handled by
            // the closure-variable reconstruction (e_int) plus the fused
            // cell pass.
            let (pl, pr) = (
                crate::eqair::eqair_pressure_at(rl, xl),
                crate::eqair::eqair_pressure_at(rr, xr),
            );
            (
                FacePrim {
                    rho: rl,
                    u: ul,
                    v: vl,
                    p: pl,
                    a: crate::eqair::eqair_sound_at(rl, xl),
                    e: xl,
                },
                FacePrim {
                    rho: rr,
                    u: ur,
                    v: vr,
                    p: pr,
                    a: crate::eqair::eqair_sound_at(rr, xr),
                    e: xr,
                },
            )
        }
    }
}

/// the face fluxes for both sweep directions, per unit area: the shared
/// machinery between the staircase update and the cut-cell update.
/// `p` and `e_int` are the precomputed cell fields of the closure.
#[allow(clippy::needless_range_loop)]
pub(crate) fn compute_axi_fluxes(
    model: ThermoModel,
    state: &ConservedState2d,
    g: &Grid2d,
    muscl: bool,
    bc: &Boundaries2d,
    p: &[f64],
    e_int: &[f64],
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(Vec<Sweep>, Vec<Sweep>), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let gamma = match model {
        ThermoModel::Perfect { gamma } => gamma,
        ThermoModel::EqAir => 1.4, // only the fallback shaping value
    };
    // the reconstructed closure variable per cell: p or e.
    let x: Vec<f64> = match model {
        ThermoModel::Perfect { .. } => p.to_vec(),
        ThermoModel::EqAir => e_int.to_vec(),
    };
    // the row and column sweeps are built by run_jobs below, sharded
    // across threads when nthreads > 1.

    // the closure variable's inflow ghosts: under equilibrium air the
    // row carries internal energy, but ghost_value answers the stored
    // pressure for supersonic inflow. invert the inflow state through
    // the closure so the ghost is the same gas the interior holds.
    let bc_x = |side: &Bc2d, model: ThermoModel| -> Option<f64> {
        match (side, model) {
            (Bc2d::SupersonicInflow { rho, p, .. }, ThermoModel::EqAir) => {
                crate::eqair::energy_from_pressure(&[*rho], &[*p])
                    .ok()
                    .and_then(|v| v.first().copied())
            }
            _ => None,
        }
    };
    let west_x = bc_x(&bc.west, model);
    let east_x = bc_x(&bc.east, model);
    let south_x_cache = bc_x(&bc.south, model);
    let north_x_cache = bc_x(&bc.north, model);

    // shard the row and column sweeps across threads when asked; each job
    // fills one sweep with no overlap, so the joined result matches the
    // serial march exactly regardless of thread count.
    let row_job = |j: usize| -> Result<Sweep, Error> {
        let (mut pr, mut pu, mut pv, mut px) =
            (vec![0.0; nx], vec![0.0; nx], vec![0.0; nx], vec![0.0; nx]);
        for i in 0..nx {
            let k = idx(i, j);
            pr[i] = state.rho[k];
            pu[i] = u[k];
            pv[i] = v[k];
            px[i] = x[k];
        }
        let y_row = g.centers_y[j * nx];
        let pad = |c: &[f64], var: usize| -> Vec<f64> {
            let mut out = vec![0.0; nx + 2];
            out[0] = west_x
                .filter(|_| var == 3)
                .unwrap_or_else(|| ghost_value(&bc.west, c[0], var, 1, y_row));
            out[nx + 1] = east_x
                .filter(|_| var == 3)
                .unwrap_or_else(|| ghost_value(&bc.east, c[nx - 1], var, 1, y_row));
            out[1..=nx].copy_from_slice(c);
            out
        };
        let flux_row = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            face_states(&padded, nx, muscl)
        };
        let (rl, rr) = flux_row(&pr, 0);
        let (ul, ur) = flux_row(&pu, 1);
        let (vl, vr) = flux_row(&pv, 2);
        let (xl, xr) = flux_row(&px, 3);
        let mut sw = Sweep::new(nx);
        for f in 0..=nx {
            let (l, r) = face_pair(
                model, rl[f], ul[f], vl[f], xl[f], rr[f], ur[f], vr[f], xr[f],
            );
            let q = hllc_flux(gamma, l, r, 0);
            sw.mass[f] = q.mass;
            sw.mx[f] = q.mx;
            sw.my[f] = q.my;
            sw.e[f] = q.e;
        }
        Ok(sw)
    };

    let col_job = |i: usize| -> Result<Sweep, Error> {
        let (mut cr, mut cu, mut cv, mut cx) =
            (vec![0.0; ny], vec![0.0; ny], vec![0.0; ny], vec![0.0; ny]);
        for j in 0..ny {
            let k = idx(i, j);
            cr[j] = state.rho[k];
            cu[j] = u[k];
            cv[j] = v[k];
            cx[j] = x[k];
        }
        let x_col = g.centers_x[i];
        let x_ghost_s = 2.0 * g.faces_x[0] - x_col;
        let x_ghost_n = 2.0 * g.faces_x[nx] - x_col;
        let pad = |c: &[f64], var: usize| -> Vec<f64> {
            let mut out = vec![0.0; ny + 2];
            out[0] = south_x_cache
                .filter(|_| var == 3)
                .unwrap_or_else(|| ghost_value(&bc.south, c[0], var, 2, x_ghost_s));
            out[ny + 1] = north_x_cache
                .filter(|_| var == 3)
                .unwrap_or_else(|| ghost_value(&bc.north, c[ny - 1], var, 2, x_ghost_n));
            out[1..=ny].copy_from_slice(c);
            out
        };
        let flux_col = |c: &[f64], var: usize| -> (Vec<f64>, Vec<f64>) {
            let padded = pad(c, var);
            face_states(&padded, ny, muscl)
        };
        let (rl, rr) = flux_col(&cr, 0);
        let (ul, ur) = flux_col(&cu, 1);
        let (vl, vr) = flux_col(&cv, 2);
        let (xl, xr) = flux_col(&cx, 3);
        let mut sw = Sweep::new(ny);
        for f in 0..=ny {
            let (l, r) = face_pair(
                model, rl[f], ul[f], vl[f], xl[f], rr[f], ur[f], vr[f], xr[f],
            );
            let q = hllc_flux(gamma, l, r, 1);
            sw.mass[f] = q.mass;
            sw.mx[f] = q.mx;
            sw.my[f] = q.my;
            sw.e[f] = q.e;
        }
        Ok(sw)
    };

    let run_jobs = |n: usize,
                    job: &(dyn Fn(usize) -> Result<Sweep, Error> + Sync)|
     -> Result<Vec<Sweep>, Error> {
        if nthreads <= 1 || n < 2 {
            return (0..n).map(job).collect();
        }
        let threads = nthreads.min(n);
        let chunk = n.div_ceil(threads);
        let mut out: Vec<Option<Sweep>> = (0..n).map(|_| None).collect();
        let mut first_err: Option<Error> = None;
        std::thread::scope(|s| {
            let mut handles = Vec::new();
            for start in (0..n).step_by(chunk) {
                let end = (start + chunk).min(n);
                let job = job;
                handles.push(s.spawn(move || (start, (start..end).map(job).collect::<Vec<_>>())));
            }
            for h in handles {
                let (start, v) = h.join().unwrap();
                for (k, sw) in v.into_iter().enumerate() {
                    match sw {
                        Ok(s) => out[start + k] = Some(s),
                        Err(e) => first_err = first_err.or(Some(e)),
                    }
                }
            }
        });
        match first_err {
            Some(e) => Err(e),
            None => Ok(out.into_iter().map(|s| s.unwrap()).collect()),
        }
    };

    let mut fx = run_jobs(ny, &row_job)?;
    let mut fy = run_jobs(nx, &col_job)?;

    // slip-wall enforcement: zero the advective flux across any face that
    // touches a solid cell. the riemann solve never runs across the wall,
    // so the interior dummy state cannot inject mass, momentum, or energy
    // into the fluid, and the fluid-side band cell keeps its projection.
    if let Some(dist) = wall {
        let (nxg, nyg) = (g.nx, g.ny);
        // the solid classification is cached per grid+body: this pass
        // reads the mask instead of re-querying the distance field, so
        // wall faces cost O(1) after the first step.
        let mask = crate::body::solid_mask(g, dist);
        let solid = |i: usize, j: usize| -> bool {
            if i >= nxg || j >= nyg {
                return false;
            }
            mask[j * nxg + i]
        };
        for j in 0..nyg {
            for f in 0..=nxg {
                // x-face f sits between cell (f-1, j) and (f, j).
                let left_solid = f > 0 && solid(f - 1, j);
                let right_solid = f < nxg && solid(f, j);
                if left_solid || right_solid {
                    fx[j].mass[f] = 0.0;
                    fx[j].mx[f] = 0.0;
                    fx[j].my[f] = 0.0;
                    fx[j].e[f] = 0.0;
                }
            }
        }
        for i in 0..nxg {
            for f in 0..=nyg {
                let below_solid = f > 0 && solid(i, f - 1);
                let above_solid = f < nyg && solid(i, f);
                if below_solid || above_solid {
                    fy[i].mass[f] = 0.0;
                    fy[i].mx[f] = 0.0;
                    fy[i].my[f] = 0.0;
                    fy[i].e[f] = 0.0;
                }
            }
        }
    }

    Ok((fx, fy))
}

/// heun two-stage step over the axisymmetric step, mirroring
/// `advance2d_rk2`'s composition over `advance2d`: the increment comes
/// from the pre-step state (the first advance call returns it).
/// the planar (non-axisymmetric) model-aware step: the same
/// closure-general face fluxes as the axi march, with the plain
/// divergence update — no annular weights, no radial pressure
/// source. this is what angled (non-zero AoA) bodies need, where
/// axisymmetry no longer holds.
pub fn advance2d_model(
    state: &mut ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let idx = |i: usize, j: usize| j * nx + i;
    let n = nx * ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
    }
    let p = thermo::pressure(model, &state.rho, &e_int, &u, &v)?;
    let a = thermo::sound_speed(model, &state.rho, &et, &u, &v)?;
    let mut smax = 0.0f64;
    for k in 0..n {
        smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
    }
    let dt = cfl * g.dx.min(g.dy) / smax.max(1e-12);

    let (fx, fy) = compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, nthreads, wall)?;

    let mut resid = 0.0f64;
    #[allow(clippy::needless_range_loop)]
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            let dxi = dt / g.dxs[i];
            let dyi = dt / g.dys[j];
            let drho = -(dxi * (fx[j].mass[i + 1] - fx[j].mass[i])
                + dyi * (fy[i].mass[j + 1] - fy[i].mass[j]));
            let dmx =
                -(dxi * (fx[j].mx[i + 1] - fx[j].mx[i]) + dyi * (fy[i].mx[j + 1] - fy[i].mx[j]));
            let dmy =
                -(dxi * (fx[j].my[i + 1] - fx[j].my[i]) + dyi * (fy[i].my[j + 1] - fy[i].my[j]));
            let de = -(dxi * (fx[j].e[i + 1] - fx[j].e[i]) + dyi * (fy[i].e[j + 1] - fy[i].e[j]));
            state.rho[k] += drho;
            state.mx[k] += dmx;
            state.my[k] += dmy;
            state.e[k] += de;
            resid = resid
                .max(drho.abs())
                .max(dmx.abs())
                .max(dmy.abs())
                .max(de.abs());
        }
    }
    Ok((dt, resid))
}

/// heun (rk2) wrapper around the planar model-aware step, matching
/// the axi march's time integration.
pub fn advance2d_model_rk2(
    state: &mut ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (dt, _) = {
        let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
        let a = thermo::sound_speed(model, &state.rho, &et, &u, &v)?;
        let mut smax = 0.0f64;
        for k in 0..g.nx * g.ny {
            smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
        }
        (cfl * g.dx.min(g.dy) / smax.max(1e-12), 0.0)
    };
    let mut s1 = state.clone();
    advance2d_model(&mut s1, g, model, cfl, muscl, bc, nthreads, wall)?;
    advance2d_model(&mut s1, g, model, cfl, muscl, bc, nthreads, wall)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt, 0.0))
}

pub fn advance2d_axi_rk2(
    state: &mut ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &Boundaries2d,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (dt, _) = {
        let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
        let n = g.nx * g.ny;
        let mut e_int = vec![0.0; n];
        for k in 0..n {
            e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        }
        let a = thermo::sound_speed(model, &state.rho, &et, &u, &v)?;
        let mut smax = 0.0f64;
        for k in 0..n {
            smax = smax.max(u[k].abs().max(v[k].abs()) + a[k]);
        }
        (cfl * g.dx.min(g.dy) / smax.max(1e-12), 0.0)
    };
    let mut s1 = state.clone();
    advance2d_axi(&mut s1, g, model, cfl, muscl, bc, nthreads, wall)?;
    advance2d_axi(&mut s1, g, model, cfl, muscl, bc, nthreads, wall)?;
    for k in 0..g.nx * g.ny {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt, 0.0))
}

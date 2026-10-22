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

#[allow(clippy::too_many_arguments)]
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
    let (_p0, _t0, mu_cells) = match model {
        ThermoModel::EqAir { .. } => crate::eqair_cea::fused_state_mu(&state.rho, &e_int, |t| {
            crate::viscous2d::sutherland_mu(mu0, t.max(1.0), t0_ref, s_param)
        }),
        ThermoModel::Park5 { .. } => {
            // species scope: the march is inviscid euler only;
            // the viscous, rans and cutcell lanes do not serve park5.
            return Err(Error::InvalidArgs);
        }
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
    // the explicit viscous stability bound is per-cell: dt from the worst
    // mu/rho in the same cell. pairing the global mu max (shock layer) with
    // the global rho min (wake floor, a different cell) overstated nu by
    // orders of magnitude and crushed dt far below the true bound.
    let nu_max = state
        .rho
        .iter()
        .zip(&mu_cells)
        .fold(0.0f64, |m, (&r, &mu)| {
            let r = r.max(1e-6);
            m.max(mu / r)
        });
    let dt_visc = 0.25 * g.dx.min(g.dy).powi(2) / nu_max.max(1e-12);
    let mut s1 = state.clone();
    let (dt1, _) =
        advance2d_model_capped(&mut s1, g, model, cfl, muscl, bc, dt_visc, nthreads, wall)?;
    crate::viscous2d::add_viscous_cells_mu(
        &mut s1,
        g,
        match model {
            ThermoModel::Perfect { gamma } => gamma,
            ThermoModel::EqAir { .. } => 1.4,
            ThermoModel::Park5 { gamma_shaping } => gamma_shaping,
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
        nthreads,
        wall,
    )?;
    let (dt2, _) =
        advance2d_model_capped(&mut s1, g, model, cfl, muscl, bc, dt_visc, nthreads, wall)?;
    crate::viscous2d::add_viscous_cells_mu(
        &mut s1,
        g,
        match model {
            ThermoModel::Perfect { gamma } => gamma,
            ThermoModel::EqAir { .. } => 1.4,
            ThermoModel::Park5 { gamma_shaping } => gamma_shaping,
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
        nthreads,
        wall,
    )?;
    for k in 0..n {
        state.rho[k] = 0.5 * state.rho[k] + 0.5 * s1.rho[k];
        state.mx[k] = 0.5 * state.mx[k] + 0.5 * s1.mx[k];
        state.my[k] = 0.5 * state.my[k] + 0.5 * s1.my[k];
        state.e[k] = 0.5 * state.e[k] + 0.5 * s1.e[k];
    }
    Ok((dt1, 0.0))
}

/// the fused (p, t, mu) pass shared by the laminar marches: per-cell
/// sutherland mu from the closure temperature (eqair: one table read).
#[allow(clippy::type_complexity)]
fn model_cell_mu(
    state: &ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    mu0: f64,
    t0_ref: f64,
    s_param: f64,
) -> Result<(Vec<f64>, Vec<f64>, Vec<f64>), Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
    }
    match model {
        ThermoModel::EqAir { .. } => {
            Ok(crate::eqair_cea::fused_state_mu(&state.rho, &e_int, |t| {
                crate::viscous2d::sutherland_mu(mu0, t.max(1.0), t0_ref, s_param)
            }))
        }
        ThermoModel::Park5 { .. } => {
            // species scope: the march is inviscid euler only;
            // the viscous, rans and cutcell lanes do not serve park5.
            Err(Error::InvalidArgs)
        }
        ThermoModel::Perfect { .. } => {
            let p = thermo::pressure(model, &state.rho, &e_int, &u, &v)?;
            let t = thermo::temperature(model, &state.rho, &p)?;
            let mu = t
                .iter()
                .map(|&tk| crate::viscous2d::sutherland_mu(mu0, tk.max(1.0), t0_ref, s_param))
                .collect::<Vec<f64>>();
            Ok((p, t, mu))
        }
    }
}

/// the thermal-diffusivity factor for the viscous dt bound: the
/// prandtl-based thermal-to-momentum diffusivity ratio, kappa ~
/// gamma/((gamma-1) pr) ~ 2.4 at pr 0.71.
fn model_kappa(pr: f64) -> f64 {
    1.4 / ((1.4 - 1.0) * pr.max(1e-6))
}

/// one local-time-stepping euler sweep of the laminar mean flow: the
/// steady-state counterpart of `advance2d_model_visc_rk2`. each cell
/// advances at its own acoustic+diffusion bound through a single
/// forward-euler stage: heun's neighbor-consistent stages are unstable
/// under per-cell increments, but at a steady state the one-stage update
/// converges to the same fixed point, and time accuracy is irrelevant
/// there. returns the smallest cell increment for progress accounting.
#[allow(clippy::too_many_arguments)]
pub fn advance2d_model_visc_lts(
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
    let (_p0, _t0, mu_cells) = model_cell_mu(state, g, model, mu0, t0_ref, s_param)?;
    let kappa_cells = vec![model_kappa(pr); n];
    let dts = model_local_dts(state, g, model, cfl, &mu_cells, &kappa_cells, wall)?;
    let (dt_min, resid) =
        advance2d_model_capped_dts(state, g, model, cfl, muscl, bc, &dts, nthreads, wall)?;
    crate::viscous2d::add_viscous_cells_mu(
        state,
        g,
        match model {
            ThermoModel::Perfect { gamma } => gamma,
            ThermoModel::EqAir { .. } => 1.4,
            ThermoModel::Park5 { gamma_shaping } => gamma_shaping,
        },
        mu0,
        &mu_cells,
        pr,
        Some(crate::viscous2d::TurbCtx {
            mu_t: &[],
            pr_t: pr,
            bc,
        }),
        crate::solver2d::TimeControl::Global(dt_min),
        Some(model),
        wall_temperature,
        nthreads,
        wall,
    )?;
    Ok((dt_min, resid))
}

#[allow(clippy::too_many_arguments)]
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
        ThermoModel::EqAir { .. } => 1.4,
        ThermoModel::Park5 { gamma_shaping } => gamma_shaping,
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
        nthreads,
        wall,
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
        nthreads,
        wall,
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
/// per-cell time increments for the model-aware (eqair or perfect) march:
/// the acoustic bound at each cell's own signal speed and width, capped by
/// an optional per-cell viscous diffusion bound. cells behind the body mask
/// step at the minimum fluid increment so no shared face over-steps.
#[allow(clippy::too_many_arguments)]
pub fn model_local_dts(
    state: &ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    cfl: f64,
    mu_cells: &[f64],
    kappa_cells: &[f64],
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<Vec<f64>, Error> {
    let (nx, ny) = (g.nx, g.ny);
    let n = nx * ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
    }
    let a = match model {
        ThermoModel::EqAir { .. } => state
            .rho
            .iter()
            .zip(&e_int)
            .map(|(&r, &e)| crate::eqair_cea::p_a_at(r, e).1)
            .collect::<Vec<f64>>(),
        ThermoModel::Perfect { gamma } => state
            .rho
            .iter()
            .zip(&e_int)
            .map(|(&r, &e)| (gamma * (gamma - 1.0) * e * r / r.max(1e-12)).sqrt())
            .collect::<Vec<f64>>(),
        ThermoModel::Park5 { .. } => {
            // the cutcell lanes do not serve park5.
            return Err(Error::InvalidArgs);
        }
    };
    let mut out = vec![0.0f64; n];
    let mut dt_min_fluid = f64::INFINITY;
    for k in 0..n {
        let (i, j) = (k % nx, k / nx);
        let h = g.dxs[i].min(g.dys[j]);
        let s = u[k].abs().max(v[k].abs()) + a[k];
        let dt_ac = cfl * h / s.max(1e-9);
        let dt = if mu_cells.is_empty() {
            dt_ac
        } else {
            let rho = state.rho[k].max(1e-6);
            let nu = (mu_cells[k] / rho) * (1.0 + kappa_cells[k]);
            let dt_dif = 0.25 * h * h / nu.max(1e-14);
            dt_ac.min(dt_dif)
        };
        let solid = wall.is_some_and(|w| w(g.centers_x[k], g.centers_y[k]) <= 0.0);
        if !solid {
            dt_min_fluid = dt_min_fluid.min(dt);
        }
        // the solid dummy state is quiescent freestream, but its faces see
        // the live flow; a solid cell stepping at the freestream acoustic
        // bound imports flux*time its fluid neighbors never agreed to.
        // the fluxes across its wall faces are zeroed, so the value only
        // matters through the face pairing with the band cells: hold it
        // at the minimum fluid increment so no shared face ever over-steps.
        let solid_dt = dt_min_fluid;
        out[k] = if solid { solid_dt } else { dt };
    }
    Ok(out)
}

/// the per-cell-dt twin of `advance2d_model_capped`: same flux sweep, each
/// cell updated at its own increment. only the converged state is
/// physically meaningful (local time stepping).
#[allow(clippy::too_many_arguments)]
fn advance2d_model_capped_dts(
    state: &mut ConservedState2d,
    g: &Grid2d,
    model: ThermoModel,
    _cfl: f64,
    muscl: bool,
    bc: &crate::solver2d::Boundaries2d,
    dts: &[f64],
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(f64, f64), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let n = nx * ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0f64; n];
    for k in 0..n {
        e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
    }
    let p = match model {
        ThermoModel::EqAir { .. } => state
            .rho
            .iter()
            .zip(&e_int)
            .map(|(&r, &e)| crate::eqair_cea::p_a_at(r, e).0)
            .collect::<Vec<f64>>(),
        ThermoModel::Perfect { gamma } => state
            .rho
            .iter()
            .zip(&e_int)
            .map(|(&r, &e)| (gamma - 1.0) * r * e)
            .collect::<Vec<f64>>(),
        ThermoModel::Park5 { .. } => {
            // the cutcell lanes do not serve park5.
            return Err(Error::InvalidArgs);
        }
    };
    let (fx, fy) =
        compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, None, nthreads, wall)?;
    // standard per-cell local stepping in residual form: each cell
    // advances at its own increment over a shared flux field. the
    // intermediate states are not conservative (a shared face moves
    // different flux*time out of each side), which is accepted for a
    // steady-state march: the converged state is the residual-zero
    // field, where the dt weighting drops out entirely. a face-paired
    // scheme (each face at the slower neighbor's increment) is exactly
    // conservative per step but its fixed point is not residual-zero,
    // so it converges to a wrong field instead: the pairing weights
    // the residual by min(dt_left, dt_right), which never vanishes
    // at a true steady state. that misformulation destroyed the first
    // two gate attempts; this is the revert to the standard form.
    let mut d_r = vec![0.0f64; n];
    let mut d_m = vec![0.0f64; n];
    let mut d_n = vec![0.0f64; n];
    let mut d_e = vec![0.0f64; n];
    let mut resids = vec![0.0f64; n];
    {
        let fxr = &fx[..];
        let fyr = &fy[..];
        let dtw = dts;
        crate::solver2d::shard_slice(&mut d_r, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let d = dtw[k] / g.dxs[i];
                let dy_ = dtw[k] / g.dys[j];
                *w = -d * (fxr[j].mass[i + 1] - fxr[j].mass[i])
                    - dy_ * (fyr[i].mass[j + 1] - fyr[i].mass[j]);
            }
        });
        crate::solver2d::shard_slice(&mut d_m, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let d = dtw[k] / g.dxs[i];
                let dy_ = dtw[k] / g.dys[j];
                *w = -d * (fxr[j].mx[i + 1] - fxr[j].mx[i])
                    - dy_ * (fyr[i].mx[j + 1] - fyr[i].mx[j]);
            }
        });
        crate::solver2d::shard_slice(&mut d_n, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let d = dtw[k] / g.dxs[i];
                let dy_ = dtw[k] / g.dys[j];
                *w = -d * (fxr[j].my[i + 1] - fxr[j].my[i])
                    - dy_ * (fyr[i].my[j + 1] - fyr[i].my[j]);
            }
        });
        crate::solver2d::shard_slice(&mut d_e, nthreads, |s, win| {
            for (w, k) in win.iter_mut().zip(s..) {
                let (i, j) = (k % nx, k / nx);
                let d = dtw[k] / g.dxs[i];
                let dy_ = dtw[k] / g.dys[j];
                *w = -d * (fxr[j].e[i + 1] - fxr[j].e[i]) - dy_ * (fyr[i].e[j + 1] - fyr[i].e[j]);
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
        for k in 0..n {
            state.rho[k] += dr[k];
            state.mx[k] += dm[k];
            state.my[k] += dn[k];
            state.e[k] += de[k];
            // local stepping lets wake cells transiently undershoot into
            // non-physical states; floor the positivity invariants the
            // closure reads need (density and internal energy), which
            // never triggers on resolved flow.
            if !(state.rho[k].is_finite() && state.rho[k] > 0.0) {
                state.rho[k] = 1e-10;
                state.mx[k] = 0.0;
                state.my[k] = 0.0;
                state.e[k] = 1e-6;
            }
            let ke = 0.5 * (state.mx[k] * state.mx[k] + state.my[k] * state.my[k])
                / state.rho[k].max(1e-12);
            if !state.e[k].is_finite() || state.e[k] - ke < 1e-9 {
                state.e[k] = ke + 1e-9;
            }
        }
    }
    let resid = resids.iter().cloned().fold(0.0f64, f64::max);
    let dt_min = dts.iter().cloned().fold(f64::INFINITY, f64::min);
    Ok((dt_min, resid))
}

#[allow(clippy::too_many_arguments)]
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
    if matches!(model, ThermoModel::Park5 { .. }) {
        // the capped-dts/cutcell lanes do not serve park5.
        return Err(Error::InvalidArgs);
    }
    // the solid interior restores to this stage-entry state after the
    // update (see the restore block at the end of the step).
    let entry_rho: Vec<f64> = state.rho.clone();
    let entry_e: Vec<f64> = state.e.clone();
    // the cached solid mask shares the smax exclusion and the restore.
    let solid_mask: Option<Vec<bool>> = wall.map(|dist| {
        let band = 1.5 * g.dx.min(g.dy);
        crate::body::body_cache(g, dist, band).mask.clone()
    });
    let (nx, ny) = (g.nx, g.ny);
    let n = nx * ny;
    let (u, v, et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let mut e_int = vec![0.0; n];
    #[allow(unused_variables)]
    let idx = |i: usize, j: usize| j * nx + i;
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
        let e_int = &e_int;
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
                let e_int = &e_int;
                sc.spawn(move || {
                    for ((w, wa), k) in pw.iter_mut().zip(aw.iter_mut()).zip(s..) {
                        let (pk, ak) = match model {
                            ThermoModel::EqAir { .. } => crate::eqair_cea::p_a_at(rho[k], e_int[k]),
                            ThermoModel::Perfect { gamma } => {
                                let pk = (gamma - 1.0) * rho[k] * e_int[k];
                                let ak = (gamma * pk / rho[k].max(1e-12)).sqrt();
                                (pk, ak)
                            }
                            ThermoModel::Park5 { .. } => {
                                unreachable!("park5 rejected at the head of the cutcell helper")
                            }
                        };
                        *w = pk;
                        *wa = ak;
                    }
                });
            }
        });
    }
    let smax = crate::solver2d::shard_rows_max(n, nthreads, |k| {
        // solid cells step at the mirror-wall pushback, which ratchets
        // their dummy velocity without any mass flux to absorb it; they
        // are restored below and their value is a boundary artifact, so
        // they cannot own the acoustic bound.
        let solid = solid_mask.as_ref().is_some_and(|m| m[k]);
        if solid {
            0.0
        } else {
            u[k].abs().max(v[k].abs()) + a[k]
        }
    });
    let dt = (cfl * g.dx.min(g.dy) / smax.max(1e-12)).min(dt_cap);

    let (fx, fy) =
        compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, None, nthreads, wall)?;

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
            // positivity floor: closed body cavities initialized with the
            // freestream pump mass cell-to-cell through their interior
            // faces (the wall faces return no mass), and an interior cell
            // can drain to the hllc floor while its neighbors stay finite.
            // an unbounded rho/u ratio then detonates the sound speed and
            // the next stage NaNs. hold drained cells at a small positive
            // density, scaling momentum/energy down with it so the
            // primitive state stays finite and the cavity can settle.
            if !(state.rho[k].is_finite() && state.rho[k] >= 1e-8) {
                // a drained or non-finite cell keeps only a vacuum seed:
                // scaling by 1e-8/rho overflows to inf (then NaN) once rho
                // underflows, so anything past the floor starts fresh.
                state.rho[k] = 1e-8;
                state.mx[k] = 0.0;
                state.my[k] = 0.0;
                state.e[k] = 0.0;
            }
        }
    }
    // restore the solid interior after every sub-stage: the wall pass
    // seals their faces but the solid-solid mirror pairs still ratchet
    // their dummy velocity stage over stage, and left alone that
    // unphysical momentum leaks into the band through the shared face
    // pairing. their state is boundary scaffolding, not flow.
    if let Some(dist) = wall {
        let band = 1.5 * g.dx.min(g.dy);
        let wall_cache = crate::body::body_cache(g, dist, band);
        for k in &wall_cache.solid_cells {
            state.rho[*k] = entry_rho[*k];
            state.mx[*k] = 0.0;
            state.my[*k] = 0.0;
            state.e[*k] = entry_e[*k];
        }
        // enclosed pockets hold the same inert reference: sealed on all
        // four faces, their fluid is a stair quantization artifact.
        for k in &wall_cache.pocket_cells {
            state.rho[*k] = entry_rho[*k];
            state.mx[*k] = 0.0;
            state.my[*k] = 0.0;
            state.e[*k] = entry_e[*k];
        }
    }
    let resid = resids.iter().cloned().fold(0.0f64, f64::max);
    Ok((dt, resid))
}

#[allow(clippy::too_many_arguments)]
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
    #[allow(unused_variables)]
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

    let (fx, fy) =
        compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, None, nthreads, wall)?;

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
    yl: &[f64; 5],
    yr: &[f64; 5],
) -> (FacePrim, FacePrim) {
    match model {
        ThermoModel::Perfect { gamma } => (
            FacePrim::perfect(rl, ul, vl, xl, gamma),
            FacePrim::perfect(rr, ur, vr, xr, gamma),
        ),
        ThermoModel::EqAir { .. } => {
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
        ThermoModel::Park5 { .. } => {
            // the closure variable `x` is the internal energy per mass and the
            // per-side composition arrives with the face from the caller's
            // species reconstruction. pressure, sound speed and specific
            // internal energy all come from the mixture closure, so a frozen
            // uniform field runs as an exact perfect-gas case.
            let face = |r: f64, u: f64, v: f64, x: f64, y: &[f64; 5]| -> FacePrim {
                FacePrim {
                    rho: r,
                    u,
                    v,
                    p: crate::species_closure::pressure_at(r, x, y),
                    a: crate::species_closure::sound_speed_at(r, x, y),
                    e: x,
                }
            };
            (face(rl, ul, vl, xl, yl), face(rr, ur, vr, xr, yr))
        }
    }
}

#[allow(clippy::too_many_arguments)]
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
    species: Option<&crate::species::SpeciesState>,
    nthreads: usize,
    wall: Option<&dyn Fn(f64, f64) -> f64>,
) -> Result<(Vec<Sweep>, Vec<Sweep>), Error> {
    let (nx, ny) = (g.nx, g.ny);
    #[allow(unused_variables)]
    let idx = |i: usize, j: usize| j * nx + i;
    let (u, v, _et) = cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let gamma = match model {
        ThermoModel::Perfect { gamma } => gamma,
        ThermoModel::EqAir { .. } => 1.4, // only the fallback shaping value
        ThermoModel::Park5 { gamma_shaping } => gamma_shaping,
    };
    // the reconstructed closure variable per cell: p or e. under park5 the
    // closure variable is the internal energy per mass, like equilibrium
    // air: the mixture closure recovers the face pressure from it.
    let x: Vec<f64> = match model {
        ThermoModel::Perfect { .. } => p.to_vec(),
        ThermoModel::EqAir { .. } | ThermoModel::Park5 { .. } => e_int.to_vec(),
    };
    // the row and column sweeps are built by run_jobs below, sharded
    // across threads when nthreads > 1.

    // the closure variable's inflow ghosts: under equilibrium air the
    // row carries internal energy, but ghost_value answers the stored
    // pressure for supersonic inflow. invert the inflow state through
    // the closure so the ghost is the same gas the interior holds.
    let bc_x = |side: &Bc2d, model: ThermoModel| -> Option<f64> {
        match (side, model) {
            (Bc2d::SupersonicInflow { rho, p, .. }, ThermoModel::EqAir { .. }) => {
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
                .unwrap_or_else(|| ghost_value(&bc.west, c[0], var, 1, y_row, px[0]));
            out[nx + 1] = east_x
                .filter(|_| var == 3)
                .unwrap_or_else(|| ghost_value(&bc.east, c[nx - 1], var, 1, y_row, px[nx - 1]));
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
        // park5: reconstruct each species' mass fraction the same way as
        // the flow variables (ghosts copy the edge cell, zero-order), then
        // stack the five species into per-face side compositions.
        let (yl_rows, yr_rows) = match species {
            Some(sp) => {
                let mut l_rows = Vec::with_capacity(crate::park::NS);
                let mut r_rows = Vec::with_capacity(crate::park::NS);
                for s in 0..crate::park::NS {
                    let mut row = vec![0.0f64; nx];
                    for i in 0..nx {
                        row[i] = sp.y[s][idx(i, j)];
                    }
                    let (l, r) = flux_row(&row, 3);
                    l_rows.push(l);
                    r_rows.push(r);
                }
                (l_rows, r_rows)
            }
            None => (Vec::new(), Vec::new()),
        };
        let y_stack = |rows: &Vec<Vec<f64>>, f: usize| -> [f64; 5] {
            let mut y = [0.0f64; 5];
            for s in 0..crate::park::NS.min(rows.len()) {
                y[s] = rows[s][f];
            }
            y
        };
        let mut sw = Sweep::new(nx);
        for f in 0..=nx {
            let (yl, yr) = if species.is_some() {
                (y_stack(&yl_rows, f), y_stack(&yr_rows, f))
            } else {
                ([0.0; 5], [0.0; 5])
            };
            let (l, r) = face_pair(
                model, rl[f], ul[f], vl[f], xl[f], rr[f], ur[f], vr[f], xr[f], &yl, &yr,
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
                .unwrap_or_else(|| ghost_value(&bc.south, c[0], var, 2, x_ghost_s, cx[0]));
            out[ny + 1] = north_x_cache.filter(|_| var == 3).unwrap_or_else(|| {
                ghost_value(&bc.north, c[ny - 1], var, 2, x_ghost_n, cx[ny - 1])
            });
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
        let (yl_cols, yr_cols) = match species {
            Some(sp) => {
                let mut l_cols = Vec::with_capacity(crate::park::NS);
                let mut r_cols = Vec::with_capacity(crate::park::NS);
                for s in 0..crate::park::NS {
                    let mut col = vec![0.0f64; ny];
                    for j in 0..ny {
                        col[j] = sp.y[s][idx(i, j)];
                    }
                    let (l, r) = flux_col(&col, 3);
                    l_cols.push(l);
                    r_cols.push(r);
                }
                (l_cols, r_cols)
            }
            None => (Vec::new(), Vec::new()),
        };
        let y_stack = |rows: &Vec<Vec<f64>>, f: usize| -> [f64; 5] {
            let mut y = [0.0f64; 5];
            for s in 0..crate::park::NS.min(rows.len()) {
                y[s] = rows[s][f];
            }
            y
        };
        let mut sw = Sweep::new(ny);
        for f in 0..=ny {
            let (yl, yr) = if species.is_some() {
                (y_stack(&yl_cols, f), y_stack(&yr_cols, f))
            } else {
                ([0.0; 5], [0.0; 5])
            };
            let (l, r) = face_pair(
                model, rl[f], ul[f], vl[f], xl[f], rr[f], ur[f], vr[f], xr[f], &yl, &yr,
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
                #[allow(clippy::redundant_locals)]
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

    // slip-wall enforcement on any advective face that touches a solid
    // cell. the interior dummy state never feeds the wall: instead the
    // fluid state is mirrored about the face and the same hllc the
    // interior faces use solves the mirrored pair. the reflected
    // riemann returns mass 0 with the stagnation pressure the wall must
    // hold, so the wall is impermeable but not force-free: pinning the
    // momentum flux to the bare cell pressure (or zeroing it) removes
    // the pushback and the stagnation band jets into the wall.
    if let Some(dist) = wall {
        let (nxg, nyg) = (g.nx, g.ny);
        // the solid classification is cached per grid+body: this pass
        // reads the mask instead of re-querying the distance field, so
        // wall faces cost O(1) after the first step. the same cache
        // carries the chord-crossed faces below.
        let band = 1.5 * g.dx.min(g.dy);
        let wall_cache = crate::body::body_cache(g, dist, band);
        let mask: &[bool] = &wall_cache.mask;
        let solid = |i: usize, j: usize| -> bool {
            if i >= nxg || j >= nyg {
                return false;
            }
            mask[j * nxg + i]
        };
        // the mirrored slip-wall riemann: reflect the fluid state about
        // the wall face and run the hllc on the mirrored pair. the
        // solve returns the wall pressure the reflection builds (always
        // at or above the cell pressure for flow approaching the wall),
        // which is the pushback a real slip wall provides. pinning the
        // momentum flux to the bare cell pressure underestimates it and
        // lets the stagnation band jet into the wall.
        let sound = |k: usize| -> f64 {
            match model {
                ThermoModel::Perfect { gamma } => {
                    (gamma * p[k].max(1e-12) / state.rho[k].max(1e-12)).sqrt()
                }
                ThermoModel::EqAir { .. } => {
                    crate::eqair_cea::sound_speed_at(state.rho[k], e_int[k])
                }
                ThermoModel::Park5 { gamma_shaping } => {
                    // the mirror-wall path only runs when the march carries a
                    // wall; park5 bodies arrive with the species march, so this arm
                    // shapes from the cell pressure like the perfect lane and
                    // the march entry rejects park5+wall before reaching here.
                    (gamma_shaping * p[k].max(1e-12) / state.rho[k].max(1e-12)).sqrt()
                }
            }
        };
        let mirror_hllc = |k: usize, axis: usize, left_fluid: bool| -> crate::hllc2d::Flux4 {
            let (rho, uu, vv, pp, aa) = (
                state.rho[k].max(1e-12),
                u[k],
                v[k],
                p[k].max(1e-12),
                sound(k),
            );
            let fluid = crate::hllc2d::FacePrim {
                rho,
                u: uu,
                v: vv,
                p: pp,
                a: aa,
                e: e_int[k],
            };
            let mirror = crate::hllc2d::FacePrim {
                rho,
                u: -uu,
                v: -vv,
                p: pp,
                a: aa,
                e: e_int[k],
            };
            let (l, r) = if left_fluid {
                (fluid, mirror)
            } else {
                (mirror, fluid)
            };
            crate::hllc2d::hllc_flux(gamma, l, r, axis)
        };
        for j in 0..nyg {
            for f in 0..=nxg {
                // x-face f sits between cell (f-1, j) and (f, j).
                let left_solid = f > 0 && solid(f - 1, j);
                let right_solid = f < nxg && solid(f, j);
                // a face the boundary crosses between two outside centers
                // is a wall too: concave wedges (udder, tail notch) leak
                // freestream straight through the surface otherwise.
                let crossed = f > 0 && f < nxg && wall_cache.wall_faces_x[j * (nxg + 1) + f];
                if left_solid || right_solid || crossed {
                    // the fluid cell with the mirror on the solid side.
                    // a solid cell flush with the domain edge leaves no
                    // fluid cell on that face: the domain-boundary flux
                    // there belongs to the outer bc, not the wall pass,
                    // so it is left untouched.
                    let side = if left_solid {
                        Some(f)
                    } else if f > 0 {
                        Some(f - 1)
                    } else {
                        None
                    };
                    if let Some(side) = side.filter(|s| *s < nxg) {
                        let q = mirror_hllc(side + j * nxg, 0, right_solid);
                        fx[j].mass[f] = q.mass;
                        fx[j].mx[f] = q.mx;
                        fx[j].my[f] = q.my;
                        fx[j].e[f] = q.e;
                    }
                }
            }
        }
        for i in 0..nxg {
            for f in 0..=nyg {
                let below_solid = f > 0 && solid(i, f - 1);
                let above_solid = f < nyg && solid(i, f);
                let crossed = f > 0 && f < nyg && wall_cache.wall_faces_y[i * (nyg + 1) + f];
                if below_solid || above_solid || crossed {
                    let side = if below_solid {
                        Some(f)
                    } else if f > 0 {
                        Some(f - 1)
                    } else {
                        None
                    };
                    if let Some(side) = side.filter(|s| *s < nyg) {
                        let q = mirror_hllc(i + side * nxg, 1, above_solid);
                        fy[i].mass[f] = q.mass;
                        fy[i].mx[f] = q.mx;
                        fy[i].my[f] = q.my;
                        fy[i].e[f] = q.e;
                    }
                }
            }
        }
    }

    Ok((fx, fy))
}

#[allow(clippy::too_many_arguments)]
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
    #[allow(unused_variables)]
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

    let (fx, fy) =
        compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, None, nthreads, wall)?;

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

#[allow(clippy::too_many_arguments)]
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

#[allow(clippy::too_many_arguments)]
/// local-time-stepping axi march: the steady-state acceleration
/// variant of `advance2d_axi`. every cell advances on its own
/// acoustic dt, so only the converged state is meaningful and the
/// reported `dt` is the field maximum (a progress marker, not a
/// physical time).
pub fn advance2d_axi_lts(
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
    let a = thermo::sound_speed(model, &state.rho, &et, &u, &v)?;
    let mut dts = vec![0.0f64; n];
    let mut dt_max = 0.0f64;
    let h = g.dx.min(g.dy);
    for (((dt_k, uk), vk), ak) in dts.iter_mut().zip(&u).zip(&v).zip(&a) {
        let s = (uk.abs().max(vk.abs()) + ak).max(1e-12);
        *dt_k = cfl * h / s;
        dt_max = dt_max.max(*dt_k);
    }
    // one forward-euler step on local dts is enough for a steady
    // march: the rk2 averaging across different local dts mixes
    // states at different pseudo-times, so the plain euler update is
    // the consistent choice here (the 2d lts lanes do the same).
    let (p, e_int) = {
        let mut e_int = vec![0.0f64; n];
        for k in 0..n {
            e_int[k] = et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k]);
        }
        (thermo::pressure(model, &state.rho, &e_int, &u, &v)?, e_int)
    };
    let (fx, fy) =
        compute_axi_fluxes(model, state, g, muscl, bc, &p, &e_int, None, nthreads, wall)?;
    let mut resid = 0.0f64;
    for j in 0..ny {
        let (rm, rp) = (g.faces_y[j], g.faces_y[j + 1]);
        let rc = 0.5 * (rm + rp);
        let dr = g.dys[j];
        for (i, fycol) in fy.iter().enumerate() {
            let k = idx(i, j);
            let dt = dts[k];
            let dxi = dt / g.dxs[i];
            let radial = |f: &Vec<f64>| (rp * f[j + 1] - rm * f[j]) / (rc * dr);
            let drho = -(dxi * (fx[j].mass[i + 1] - fx[j].mass[i]) + dt * radial(&fycol.mass));
            let dmx = -(dxi * (fx[j].mx[i + 1] - fx[j].mx[i]) + dt * radial(&fycol.mx));
            let dmy =
                -(dxi * (fx[j].my[i + 1] - fx[j].my[i]) + dt * radial(&fycol.my)) + dt * p[k] / rc;
            let de = -(dxi * (fx[j].e[i + 1] - fx[j].e[i]) + dt * radial(&fycol.e));
            // positivity guard for the local step: a fast cell that
            // outruns a strong transient can drive its density or
            // internal energy negative in one euler step; scaling the
            // increment back to keep the state physical preserves the
            // march instead of aborting it.
            let rho_new = state.rho[k] + drho;
            let e_int_new = state.e[k] + de
                - 0.5 * (state.mx[k] + dmx).powi(2) / rho_new
                - 0.5 * (state.my[k] + dmy).powi(2) / rho_new;
            if rho_new > 0.0 && e_int_new > 0.0 {
                state.rho[k] = rho_new;
                state.mx[k] += dmx;
                state.my[k] += dmy;
                state.e[k] += de;
                resid = resid
                    .max(drho.abs())
                    .max(dmx.abs())
                    .max(dmy.abs())
                    .max(de.abs());
            } else {
                // halve the increment once; if still bad, skip the
                // cell this pseudo-step.
                let rho_half = state.rho[k] + 0.5 * drho;
                let e_half = state.e[k] + 0.5 * de
                    - 0.5 * (state.mx[k] + 0.5 * dmx).powi(2) / rho_half
                    - 0.5 * (state.my[k] + 0.5 * dmy).powi(2) / rho_half;
                if rho_half > 0.0 && e_half > 0.0 {
                    state.rho[k] = rho_half;
                    state.mx[k] += 0.5 * dmx;
                    state.my[k] += 0.5 * dmy;
                    state.e[k] += 0.5 * de;
                }
            }
        }
    }
    Ok((dt_max, resid))
}

#[allow(clippy::too_many_arguments)]
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

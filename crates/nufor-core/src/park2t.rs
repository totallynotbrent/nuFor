//! park two-temperature model: a live vibrational energy field e_v (J/kg)
//! advected with the euler mass fluxes and relaxed toward the local
//! translational temperature via landau-teller (millikan-white with the
//! park high-temperature cutoff). the translational pool sees the
//! vibrational draw, so total energy is conserved cell-for-cell.
//!
//! the model is operator-split: the euler march (serial rust or fortran
//! mpi, unchanged) advances the conserved field; this module then
//! advects e_v with the same face mass fluxes and integrates the
//! relaxation ode backward-euler, which is unconditionally stable and
//! exact for the frozen-t limit. the vibrational pool feeds back
//! through the energy split only; the reaction-temperature coupling
//! (t_rxn = sqrt(t * t_v)) activates with the species march.

use crate::error::Error;
use crate::state2d::ConservedState2d;
use crate::Grid2d;

/// n2 molar mass, kg/mol.
pub const M_N2: f64 = 28.016e-3;
/// n2 characteristic vibrational temperature, K (park 1990).
pub const THETA_N2: f64 = 3393.0;
/// specific gas constant of n2, J/(kg K).
pub const R_N2: f64 = 8.314462 / M_N2;

/// sutherland viscosity for n2 (white): mu0 = 1.663e-5 Pa s at
/// t0 = 273.15 K, s = 106.7 K. used by the millikan-white partner
/// density when the cutoff term needs the collision rate.
pub fn sutherland_mu_n2(t: f64) -> f64 {
    1.663e-5 * (t / 273.15).powf(1.5) * (273.15 + 106.7) / (t + 106.7)
}

/// harmonic-oscillator vibrational energy per mass at vibrational
/// temperature t_v, J/kg.
pub fn e_vib_of_t_v(t_v: f64) -> f64 {
    let x = THETA_N2 / t_v.max(1.0);
    if x > 40.0 {
        return 0.0;
    }
    let ex = x.exp();
    R_N2 * THETA_N2 / (ex - 1.0).max(1e-12)
}

/// invert e_vib_of_t_v: vibrational temperature for a stored
/// vibrational energy. bisection over [1, 60000] K; the caloric form
/// is monotone in t_v so the bracket is exact.
pub fn t_vib_of_e_v(e_v: f64) -> f64 {
    if e_v <= 0.0 {
        return 1.0;
    }
    let (mut lo, mut hi) = (1.0f64, 60000.0f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if e_vib_of_t_v(mid) < e_v {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// millikan-white vibrational relaxation time for n2-n2 collisions
/// with park's collision-limited cutoff. t and t_v in K, p in Pa,
/// number density n in 1/m^3; returns seconds.
pub fn tau_vt_n2(t: f64, p_pa: f64, n: f64) -> f64 {
    // reduced mass of the colliding pair n2-n2, g/mol.
    let mu_gmol = M_N2 * 1e3 / 2.0;
    let inv_t = 1.0 / t.max(1.0);
    let a_s = 1.16e-3 * mu_gmol.sqrt() * THETA_N2.powf(4.0 / 3.0);
    let b_s = 0.015 * mu_gmol.powf(0.25);
    let expo = (a_s * (inv_t.cbrt() - b_s) - 18.42).clamp(-50.0, 50.0);
    let tau_mw = expo.exp() / p_pa.max(1.0);
    // park cutoff: sigma_v = 1e-21 * t^2 m^2, mean thermal speed of the
    // reduced pair.
    let r_pair = 8.314462 / (mu_gmol * 1e-3);
    let c_bar = (8.0 * r_pair * t / std::f64::consts::PI).sqrt();
    let tau_cut = 1.0 / (1e-21 * t * t * n.max(1.0) * c_bar);
    tau_mw + tau_cut
}

/// the vibrational field carried per run: one e_v per cell, J/kg.
#[derive(Clone, Debug)]
pub struct Park2tState {
    pub e_v: Vec<f64>,
}

impl Park2tState {
    pub fn from_uniform(g: &Grid2d, e_v0: f64) -> Self {
        Self {
            e_v: vec![e_v0; g.nx * g.ny],
        }
    }
}

/// one advect+relax substep on the vibrational field, drawn against the
/// post-euler state. the advection reuses the euler update's face mass
/// fluxes through a donor-cell (first-order upwind) pass; the
/// relaxation is backward-euler with the energy draw removed from the
/// translational pool so total energy is conserved exactly.
///
/// gamma feeds the pressure recovery on the translational energy.
#[allow(clippy::too_many_arguments)]
pub fn park2t_substep(
    state: &mut ConservedState2d,
    park: &mut Park2tState,
    g: &Grid2d,
    gamma: f64,
    dt: f64,
    solid: &[bool],
    muscl_e_v: bool,
) -> Result<(), Error> {
    let (nx, ny) = (g.nx, g.ny);
    let n = nx * ny;
    if park.e_v.len() != n {
        return Err(Error::InvalidArgs);
    }
    let idx = |i: usize, j: usize| j * nx + i;

    // translational pressure and temperature from the post-euler state:
    // ein here is the full internal pool (e_v included), so the split
    // removes e_v before the gas constant sees it.
    let mut p = vec![0.0f64; n];
    let mut t = vec![0.0f64; n];
    for k in 0..n {
        let rho = state.rho[k];
        let u = state.mx[k] / rho.max(1e-12);
        let v = state.my[k] / rho.max(1e-12);
        let et = state.e[k] / rho.max(1e-12);
        let ein_full = et - 0.5 * (u * u + v * v);
        let ein_tr = (ein_full - park.e_v[k]).max(1.0);
        p[k] = (gamma - 1.0) * rho * ein_tr;
        t[k] = p[k] / (rho.max(1e-12) * R_N2);
    }

    // donor-cell advection of e_v on the same face mass fluxes the
    // euler step used: face rho*u from the cell averages. solid faces
    // carry no mass flux (the euler closure zeroes them), so the
    // upwind donor there is the same cell: no transfer.
    let mut e_new = park.e_v.clone();
    for j in 0..ny {
        for i in 0..nx {
            if solid[idx(i, j)] {
                continue;
            }
            let k = idx(i, j);
            let rho = state.rho[k].max(1e-12);
            let u = state.mx[k] / rho;
            let v = state.my[k] / rho;
            let mut d = 0.0f64;
            // x faces
            let ev_here = park.e_v[k];
            if i > 0 && !solid[idx(i - 1, j)] {
                let f = if u >= 0.0 {
                    rho * u * park.e_v[idx(i - 1, j)]
                } else {
                    rho * u * ev_here
                };
                d -= f / g.dx;
            }
            if i + 1 < nx && !solid[idx(i + 1, j)] {
                let f = if u >= 0.0 {
                    rho * u * ev_here
                } else {
                    rho * u * park.e_v[idx(i + 1, j)]
                };
                d += f / g.dx;
            }
            // y faces
            if j > 0 && !solid[idx(i, j - 1)] {
                let f = if v >= 0.0 {
                    rho * v * park.e_v[idx(i, j - 1)]
                } else {
                    rho * v * ev_here
                };
                d -= f / g.dy;
            }
            if j + 1 < ny && !solid[idx(i, j + 1)] {
                let f = if v >= 0.0 {
                    rho * v * ev_here
                } else {
                    rho * v * park.e_v[idx(i, j + 1)]
                };
                d += f / g.dy;
            }
            let _ = muscl_e_v; // reserved for the second-order lane
            e_new[k] = (park.e_v[k] + dt * d).max(0.0);
        }
    }
    park.e_v = e_new;

    // backward-euler relaxation with the translational draw.
    // e_new = (e_v + (dt/tau) e_eq(t)) / (1 + dt/tau)
    // the vibrational pool lives INSIDE the conserved total (ein_tr =
    // ein_full - e_v), so the draw moves energy between pools inside
    // state.e and the total is conserved by construction.
    for k in 0..n {
        if solid[k] {
            continue;
        }
        let rho = state.rho[k].max(1e-12);
        let n_den = rho / M_N2 * 6.022e23;
        let tau = tau_vt_n2(t[k], p[k], n_den);
        let e_eq = e_vib_of_t_v(t[k]);
        let w = dt / tau.max(1e-9);
        let e_new = (park.e_v[k] + w * e_eq) / (1.0 + w);
        park.e_v[k] = e_new;
    }
    Ok(())
}

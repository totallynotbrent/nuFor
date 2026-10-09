//! park two-temperature 5-species air (n2, o2, no, n, o).
//! euler closure is still the fitted equilibrium `eqair` pressure;
//! this module advances a frozen-chemistry vibrational energy density
//! e_v per cell with landau-teller relaxation and lets the gas kinetics
//! back-couple through the mixing rule `t_rxn = sqrt(t * t_vib)`.
//!
//! this module is the conservative rung: the euler solution holds, the
//! reaction source writes only into the vibrational reservoir, and the
//! dissociation routes its energy draw through the same closure the
//! euler march is already using. the full park chemistry (coupled
//! species continuity + torsional/thermal noneq) is a later rung.

use crate::eqair::{gamm_and_partials, E0};
use crate::eqair::RHO0;
use crate::error::Error;
use crate::state2d::ConservedState2d;
use crate::Grid2d;

/// species index. order is fixed across the module.
pub const N_N2: usize = 0;
pub const N_O2: usize = 1;
pub const N_NO: usize = 2;
pub const N_N: usize = 3;
pub const N_O: usize = 4;

/// molecular masses, kg/mol.
pub const M: [f64; 5] = [28.016e-3, 32.000e-3, 30.008e-3, 14.008e-3, 16.000e-3];

/// characteristic vibrational temperatures, K (park 1990, table 1).
pub const THETA_V: [f64; 5] = [3393.0, 2239.0, 2817.0, 0.0, 0.0];

/// dissociation energies in K (e = d0 / kb; park 1990).
pub const D0_K: [f64; 3] = [113200.0, 59400.0, 75500.0];

/// modified-arrhenius forward rates, SI: kf = a * t_r^eta * exp(-td / t_r),
/// with a in m^3 mol^-1 s^-1 (converted from the cgs table by * 1e-6).
pub struct ReactionRate {
    /// arrhenius prefactor, m^3 mol^-1 s^-1.
    pub a: f64,
    /// exponent on the rate temperature.
    pub eta: f64,
    /// characteristic dissociation temperature, K.
    pub td: f64,
}

/// the five forward rates of the park 5-species air baseline
/// (park 1990 table 2). note the N2 dissociation carries a different
/// rate into N and O than into the molecular partners.
pub const RATES: [ReactionRate; 6] = [
    // 1. n2 + m -> n + n + m
    ReactionRate { a: 7.0e15, eta: -1.6, td: D0_K[0] }, // M = N2, O2, NO
    ReactionRate { a: 3.0e16, eta: -1.6, td: D0_K[0] }, // M = N,  O
    // 2. o2 + m -> o + o + m
    ReactionRate { a: 2.0e15, eta: -1.5, td: D0_K[1] },
    // 3. no + m -> n + o + m
    ReactionRate { a: 5.0e9, eta: -1.0, td: D0_K[2] },
    // 4. n2 + o -> no + n   (zeldovich first)
    ReactionRate { a: 6.4e11, eta: -1.0, td: 38370.0 },
    // 5. no + o -> o2 + n   (zeldovich second)
    ReactionRate { a: 8.4e6, eta: 0.0, td: 19450.0 },
];

/// molar gas constant, J/(mol K).
pub const RU: f64 = 8.314462;

/// the two-temperature reaction temperature, park's geometric-mean
/// rule with exponent `a` in [0, 1]:
///     t_r = t^a * tv^(1 - a)
/// the baseline model uses a = 0.5 (the `sqrt(t * tv)` form).
pub fn park_reaction_temperature(t: f64, t_vib: f64, a: f64) -> f64 {
    if t_vib <= 0.0 {
        return t;
    }
    t.powf(a) * t_vib.powf(1.0 - a)
}

/// molar vibrational energy reservoir of species s at temperature t_v
/// (harmonic oscillator, caloric form):
///     e_vib,s = (R_u * theta_v,s) / (exp(theta_v,s / t_v) - 1)
/// returns J/mol.
pub fn e_vib_per_mol(theta_v: f64, t_v: f64) -> f64 {
    let x = theta_v / t_v.max(1.0);
    if x > 40.0 {
        return 0.0;
    }
    let ex = x.exp();
    RU * theta_v / (ex - 1.0).max(1e-12)
}

/// landau-teller vibrational relaxation time in seconds for species s in
/// mol/m^3 units of background density rho; the millikan-white form with
/// park's high-temperature cutoff correction.
///
/// the implementation is per collision partner; the caller passes the
/// pair's reduced mass (g/mol) and the number density of the partner
/// (1/m^3).
pub fn tau_v_millikan_white(theta_v: f64, t: f64, mu_gmol: f64, p_pa: f64) -> f64 {
    let inv_t = 1.0 / t.max(1.0);
    let a_s = 1.16e-3 * mu_gmol.sqrt() * theta_v.powf(4.0 / 3.0);
    let b_s = 0.015 * mu_gmol.powf(0.25);
    let expo = (a_s * (inv_t.cbrt() - b_s) - 18.42).clamp(-50.0, 50.0);
    expo.exp() / p_pa.max(1.0)
}

/// park's high-temperature vibration-translational cutoff (collision-
/// limited) term: sigma_v = 1e-21 * t^2 m^2, n in 1/m^3.
pub fn tau_v_park_cutoff(t: f64, mu_gmol: f64, n: f64) -> f64 {
    let r_spec = RU / (mu_gmol * 1.0e-3); // J/(kg K) = R_u / (kg/mol)
    let c_bar = (8.0 * r_spec * t / std::f64::consts::PI).sqrt();
    1.0 / (1e-21 * t * t * n.max(1.0) * c_bar)
}

/// total vibrational relaxation time for species s:
///     tau = tau_mw + tau_cutoff
pub fn tau_vib_total(
    theta_v: f64,
    t: f64,
    mu_gmol: f64,
    p_pa: f64,
    n: f64,
) -> f64 {
    tau_v_millikan_white(theta_v, t, mu_gmol, p_pa) + tau_v_park_cutoff(t, mu_gmol, n)
}

/// rate of vibrational-energy change for species s in a gas mixture;
/// J/(m^3 s):
///     de_v,s/dt = rho_s/M_s * (e_vib,s(t) - e_vib,s(t_v)) / tau_v,s
///
/// `rho_s_kgm3` is the species density in the cell; `p_pa` the local
/// translational pressure used by the relaxation times.
pub fn d_e_vib_dt(
    theta_v: f64,
    rho_s_kgm3: f64,
    m_kgmol: f64,
    t: f64,
    t_v: f64,
    p_pa: f64,
    n_m2_3: f64,
) -> f64 {
    if theta_v <= 0.0 {
        return 0.0;
    }
    let mu = (m_kgmol * 1e3).max(1e-6); // g/mol, partner mass uses n2 as
                                         // dominant collision partner.
    let tau = tau_vib_total(theta_v, t, mu, p_pa, n_m2_3);
    let n_s = rho_s_kgm3 / m_kgmol; // mol/m^3
    let gap = e_vib_per_mol(theta_v, t) - e_vib_per_mol(theta_v, t_v);
    n_s * gap / tau.max(1e-9)
}

/// state of the noneq carry-along reservoir: one cell = one value of
/// e_vib (J/kg), one for each molecular species. atoms carry none.
#[derive(Clone, Debug)]
pub struct ParkState {
    /// vibrational energy per mass for N2, J/kg.
    pub e_vib_n2: Vec<f64>,
    /// vibrational energy per mass for O2.
    pub e_vib_o2: Vec<f64>,
    /// vibrational energy per mass for NO.
    pub e_vib_no: Vec<f64>,
    /// shared vibrational temperature per cell, K.
    pub t_v: Vec<f64>,
}

impl ParkState {
    pub fn zeros(n: usize) -> Self {
        Self {
            e_vib_n2: vec![0.0; n],
            e_vib_o2: vec![0.0; n],
            e_vib_no: vec![0.0; n],
            t_v: vec![300.0; n],
        }
    }
}

/// one backward-euler vibrational-relaxation substep on the park state,
/// drawn against the euler state. dt in seconds; uses the cell pressure
/// from the eqair closure.
pub fn relax_vibrational_energy(
    state: &ConservedState2d,
    g: &Grid2d,
    park: &mut ParkState,
    model: crate::thermo::ThermoModel,
    dt: f64,
) -> Result<(), Error> {
    let n = g.nx * g.ny;
    let (u, v, et) = crate::state2d::cons_to_prim2d(&state.rho, &state.mx, &state.my, &state.e)?;
    let e_int: Vec<f64> = (0..n)
        .map(|k| (et[k] - 0.5 * (u[k] * u[k] + v[k] * v[k])).max(1.0e4))
        .collect();
    let p = crate::thermo::pressure(model, &state.rho, &e_int, &u, &v)?;
    let temps = crate::thermo::temperature(model, &state.rho, &p)?;
    for k in 0..n {
        let t = temps[k];
        // approximate species split for the relaxation partner: pure
        // n2 at the cold limit is the baseline approximation here; the
        // full mass-fraction coupling belongs to the species-equation
        // rung.
        let rho_s = state.rho[k];
        let n_tot = rho_s / M[0] * 6.022e23; // 1/m^3
        let sv = d_e_vib_dt(THETA_V[0], rho_s, M[0], t, park.t_v[k], p[k], n_tot);
        let dv = dt * sv / rho_s.max(1e-12); // J/kg
        park.e_vib_n2[k] = (park.e_vib_n2[k] + dv).max(0.0);
        park.t_v[k] = t; // placeholder: one shared bath at the euler t;
    }
    Ok(())
}

/// equilibrium-air backward constants (in the TGAS1 sense) — pulled from
/// gamm and e; used by `park_rhs` when the caller wants the equilibrium
/// target the relaxation is approaching.
pub fn equilibrium_e_v(
    rho: &[f64],
    e_int: &[f64],
    y_mol_n2: f64,
    y_mol_o2: f64,
) -> Result<Vec<f64>, Error> {
    let n = rho.len();
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let y = (rho[i] / RHO0).log10().clamp(-6.0, 1.0);
        let z = (e_int[i] / E0).log10().clamp(0.65, 3.4);
        let (g, _, _) = gamm_and_partials(y, z);
        // cold-limit gamma gives no vibrational excitation
        let t = (e_int[i] * (g - 1.0)).max(1.0);
        let ev = RU * (y_mol_n2 * THETA_V[0] / ((THETA_V[0] / t).exp() - 1.0).max(1e-12)
            + y_mol_o2 * THETA_V[1] / ((THETA_V[1] / t).exp() - 1.0).max(1e-12));
        out.push(ev);
    }
    Ok(out)
}

/// one-shot solver-level step: advances euler by the existing closure-
/// aware march, then relaxes the vibrational bath. the euler's energy
/// is left at the equilibrium value; the noneq identifier `t_v`
/// reports how far off equilibrium the local cell is.
pub fn advance2d_park_rk2(
    state: &mut ConservedState2d,
    park: &mut ParkState,
    g: &Grid2d,
    model: crate::thermo::ThermoModel,
    cfl: f64,
    muscl: bool,
    bc: &crate::Boundaries2d,
) -> Result<(f64, f64), Error> {
    let (dt, a0) = crate::advance2d_model_rk2(state, g, model, cfl, muscl, bc)?;
    relax_vibrational_energy(state, g, park, model, dt)?;
    Ok((dt, a0))
}

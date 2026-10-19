//! park 5-species mixture closure (rung 8): pressure, temperature,
//! and sound speed from (rho, e_internal, Y) with species-dependent
//! gas constants and formation energies.
//!
//! one-temperature thermal model, vibrational energy frozen this
//! rung (it rides with the park2t machinery when two-temperature
//! lands). per-species thermal internal energy per mass:
//!
//! - atoms (N, O): e_th = (3/2) R_i T
//! - molecules (N2, O2, NO): e_th = (5/2) R_i T
//!
//! and the mixture's internal energy per mass carries the chemical
//! part through formation enthalpies referenced to the atom ground
//! states: h_f(N2) = -D0(N2), h_f(O2) = -D0(O2), h_f(NO) = -D0(NO),
//! h_f(N) = h_f(O) = 0, all per unit mass, converted from the D0
//! characteristic temperatures through D0 * kb / m_species.
//!
//! solving e(rho, T, Y) for T is a linear solve per cell: e = a*T + b
//! with a = sum Y_i cv_i and b = sum Y_i h_f_i, so T = (e - b)/a and
//! the pressure follows as p = rho R_mix T.

use crate::park::{D0_K, M, NS, RU};

/// per-species gas constant, J/(kg K).
pub const R_SPECIES: [f64; NS] = [RU / M[0], RU / M[1], RU / M[2], RU / M[3], RU / M[4]];

/// per-species constant-v (thermal only, vibration excluded),
/// J/(kg K): molecules 5/2 R, atoms 3/2 R.
pub const CV_SPECIES: [f64; NS] = [
    2.5 * RU / M[0],
    2.5 * RU / M[1],
    2.5 * RU / M[2],
    1.5 * RU / M[3],
    1.5 * RU / M[4],
];

/// per-species formation enthalpy per mass, J/kg, on the 0-K
/// atom-positive convention of TP-2867 Table I: the atoms carry the
/// dissociation energy (N: 33.6 MJ/kg, O: 15.4), NO carries the Hess
/// cycle residual, and the ground-state molecules are zero.
pub const HF_SPECIES: [f64; NS] = [
    0.0,
    0.0,
    RU * (0.5 * D0_K[0] + 0.5 * D0_K[1] - D0_K[2]) / M[2],
    0.5 * RU * D0_K[0] / M[3],
    0.5 * RU * D0_K[1] / M[4],
];

/// mixture thermal coefficients for one cell: (a, b) with
/// e = a*T + b. `a` is the mixture cv (thermal), `b` the formation
/// energy per mass.
fn mixture_coeffs(y_cell: &[f64]) -> (f64, f64) {
    let mut a = 0.0f64;
    let mut b = 0.0f64;
    for i in 0..NS {
        a += y_cell[i] * CV_SPECIES[i];
        b += y_cell[i] * HF_SPECIES[i];
    }
    (a, b)
}

/// mixture gas constant per cell, J/(kg K).
pub fn mix_gas_constant(y_cell: &[f64]) -> f64 {
    let mut r = 0.0f64;
    for i in 0..NS {
        r += y_cell[i] * R_SPECIES[i];
    }
    r
}

/// temperature from (e_internal, Y) for one cell, K. thermal energy
/// only: the formation part subtracts before the linear solve.
pub fn temperature_at(e_int: f64, y_cell: &[f64]) -> f64 {
    let (a, b) = mixture_coeffs(y_cell);
    // cold mixtures can transiently carry e below the formation
    // floor during a shock formation step; floor the thermal part
    // rather than fail the step, the same policy the eqair clamp
    // uses for out-of-window cells.
    ((e_int - b) / a).max(50.0)
}

/// internal energy per mass at (T, Y) for one cell: the forward
/// direction of `temperature_at`, the mixture's thermal-plus-
/// formation energy the conserved state carries.
pub fn internal_energy_at(t: f64, y_cell: &[f64]) -> f64 {
    let (a, b) = mixture_coeffs(y_cell);
    a * t + b
}

/// pressure from (rho, e_internal, Y) for one cell, Pa.
pub fn pressure_at(rho: f64, e_int: f64, y_cell: &[f64]) -> f64 {
    let t = temperature_at(e_int, y_cell);
    rho * mix_gas_constant(y_cell) * t
}

/// frozen-composition sound speed: gamma_mix from cp/cv with the
/// thermal cp = cv + R per species (vibration excluded), the
/// standard one-temperature frozen gamma for the 5-species set.
pub fn sound_speed_at(rho: f64, e_int: f64, y_cell: &[f64]) -> f64 {
    let mut cv = 0.0f64;
    let mut cp = 0.0f64;
    for i in 0..NS {
        cv += y_cell[i] * CV_SPECIES[i];
        cp += y_cell[i] * (CV_SPECIES[i] + R_SPECIES[i]);
    }
    let gamma = (cp / cv).max(1.01);
    let p = pressure_at(rho, e_int, y_cell).max(1e-9);
    (gamma * p / rho).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn air_y() -> [f64; NS] {
        [0.767, 0.233, 0.0, 0.0, 0.0]
    }

    /// cold air must reproduce the perfect-gas values for the
    /// 5-species subset (no argon): R = 288.2 J/(kg K), cv =
    /// 2.5 R = 720.4, gamma = 1.4, p = rho R T.
    #[test]
    fn cold_air_matches_perfect_gas() {
        let y = air_y();
        let r = mix_gas_constant(&y);
        assert!((r - 288.17).abs() < 0.05, "R = {r}");
        let (a, b) = mixture_coeffs(&y);
        assert!((a - 720.4).abs() < 0.5, "cv = {a}");
        // cold air is all ground-state molecules: the atom-positive
        // formation offset is exactly zero, so e = a*T is the sensible
        // energy and the perfect-gas correspondence is exact.
        assert!(b.abs() < 1e-9, "b = {b}");
        let t = 300.0;
        let e = a * t + b;
        let t_back = temperature_at(e, &y);
        assert!((t_back - t).abs() < 1e-9);
        let p = pressure_at(1.2, e, &y);
        let p_pg = 1.2 * r * t;
        assert!((p - p_pg).abs() / p_pg < 1e-12, "p = {p}");
        let a_snd = sound_speed_at(1.2, e, &y);
        let gamma = a_snd * a_snd / (r * t);
        assert!((gamma - 1.4).abs() < 1e-9, "gamma = {gamma}");
    }

    /// the pure-atom limit: gamma must land at 5/3 and R at RU/14 for N.
    #[test]
    fn atomic_limit() {
        let mut y = [0.0; NS];
        y[3] = 1.0; // pure N
        let r = mix_gas_constant(&y);
        assert!((r - RU / M[3]).abs() < 1e-12);
        let (a, _b) = mixture_coeffs(&y);
        assert!((a - 1.5 * RU / M[3]).abs() < 1e-12);
        // the closure energy includes the atom's formation enthalpy, so
        // build it through internal_energy_at, not the bare thermal part.
        let e = internal_energy_at(9000.0, &y);
        let a_snd = sound_speed_at(1e-4, e, &y);
        let gamma = a_snd * a_snd / (RU / M[3] * 9000.0);
        assert!((gamma - 5.0 / 3.0).abs() < 1e-9, "gamma = {gamma}");
    }

    /// dissociation absorbs sensible energy: at fixed e and rho,
    /// partly-dissociated nitrogen runs much cooler than pure N2.
    #[test]
    fn dissociation_cools_the_mixture() {
        let mut y2 = [0.0; NS];
        y2[0] = 1.0; // pure N2
        let (a2, b2) = mixture_coeffs(&y2);
        let mut y1 = [0.0; NS];
        y1[0] = 0.5;
        y1[3] = 0.5; // half-dissociated
        let (_a1, b1) = mixture_coeffs(&y1);
        let e = b2 + a2 * 12000.0; // pure N2 at 12000 K thermal
        let t2 = temperature_at(e, &y2);
        let t1 = temperature_at(e, &y1);
        assert!((t2 - 12000.0).abs() < 1e-6);
        assert!(t1 < 7000.0, "t1 = {t1}");
        // and the formation floor moved: b1 > b2 (less bound energy).
        assert!(b1 > b2);
    }
}

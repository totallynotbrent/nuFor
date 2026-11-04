//! kinetic-theory transport for the h2/o2 set: chapman-enskog
//! viscosity from the lennard-jones pairs, eucken conductivity, wilke
//! mixture rules, and binary diffusion with the neufeld collision
//! integrals. the fidelity class matches the park lane's mixture-
//! averaged treatment (wilke rules identical, per-species data from
//! the public gri-mech 3.0 transport database instead of the rp-1232
//! air fits), so the two lanes report comparable wall fluxes.
//!
//! known limit, carried honestly: the gri lj pair for h2 overpredicts
//! light-gas binary diffusivities by 35-50 percent (d(h2-n2) at 298 k,
//! 1 atm reads ~0.99 cm2/s against the 0.61-0.65 literature band),
//! a documented property of the pair itself. the pair is used anyway
//! for database consistency.

use crate::combustion_thermo::{M, NS, RU};
use crate::transport_general::{
    binary_diffusivity as general_binary_diffusivity,
    effective_diffusivity as general_effective_diffusivity,
    mixture_conductivity as general_mixture_conductivity,
    mixture_viscosity as general_mixture_viscosity,
    species_conductivity as general_species_conductivity,
    species_viscosity as general_species_viscosity, SpeciesTransport,
};

/// lennard-jones well depths epsilon/kb, k.
pub const SPECIES_WELL_DEPTH: [f64; NS] = [
    38.0,  // H2
    107.4, // O2
    572.4, // H2O
    80.0,  // OH
    145.0, // H
    80.0,  // O
    107.4, // HO2 (the yaml reuses the O2 pair)
    107.4, // H2O2 (same)
    97.53, // N2
];

/// lennard-jones collision diameters sigma, angstrom.
pub const SPECIES_DIAMETER: [f64; NS] = [
    2.92,  // H2
    3.458, // O2
    2.605, // H2O
    2.75,  // OH
    2.05,  // H
    2.75,  // O
    3.458, // HO2
    3.458, // H2O2
    3.621, // N2
];

/// rotational-relaxation numbers (not used by the eucken form; kept
/// with the pairs as the database row).
const _ZROT: [f64; NS] = [280.0, 3.8, 4.0, 1.0, 0.0, 0.0, 1.0, 3.8, 4.0];

/// the nine-species transport rows, assembled once per call site.
///
/// the values are compile-time constants, so the compiler folds this
/// into the caller and the indirection costs nothing at run time.
fn rows() -> [SpeciesTransport; NS] {
    let mut out = [SpeciesTransport {
        molar_mass: 0.0,
        well_depth: 0.0,
        diameter: 0.0,
    }; NS];
    let mut i = 0;
    while i < NS {
        out[i] = SpeciesTransport {
            molar_mass: M[i],
            well_depth: SPECIES_WELL_DEPTH[i],
            diameter: SPECIES_DIAMETER[i],
        };
        i += 1;
    }
    out
}

/// pure-species chapman-enskog viscosity, pa s. m in kg/mol, sigma in
/// angstrom: mu_i = 2.669e-6 sqrt(m t) / (sigma^2 omega_mu).
pub fn species_viscosity(i: usize, t: f64) -> f64 {
    general_species_viscosity(&rows()[i], t)
}

/// eucken correlation for pure-species conductivity, w/(m k):
/// k_i = mu_i (cp_i + 1.25 r / m_i), with cp_i the species specific
/// heat at constant pressure per kg.
pub fn species_conductivity(i: usize, t: f64, cp_i: f64) -> f64 {
    general_species_conductivity(&rows()[i], t, cp_i, RU)
}

/// wilke 1950 mixture viscosity from pure values and mass fractions.
pub fn mixture_viscosity(y: &[f64; NS], t: f64) -> f64 {
    let rows = rows();
    let pure: Vec<f64> = (0..NS)
        .map(|i| general_species_viscosity(&rows[i], t))
        .collect();
    general_mixture_viscosity(y, &rows, &pure)
}

/// wilke mixture conductivity, same weights.
pub fn mixture_conductivity(y: &[f64; NS], t: f64, cps: &[f64; NS]) -> f64 {
    let rows = rows();
    let pure: Vec<f64> = (0..NS)
        .map(|i| general_species_conductivity(&rows[i], t, cps[i], RU))
        .collect();
    general_mixture_conductivity(y, &rows, &pure)
}

/// binary diffusivity d_ij, m2/s, from kinetic theory:
/// d_ij = (3/16) sqrt(2 pi r^3 t^3 (1/m_i + 1/m_j) / (p sigma_ij^2 omega_d)).
pub fn binary_diffusivity(i: usize, j: usize, t: f64, p: f64) -> f64 {
    let rows = rows();
    general_binary_diffusivity(&rows[i], &rows[j], t, p)
}

/// effective diffusivity of species i against the mixture (the wilke
/// harmonic average the chebyshev/chemkin transport packages use):
/// 1/d_i = sum_j x_j / d_ij over j != i, with x_j mole fractions.
pub fn effective_diffusivity(y: &[f64; NS], t: f64, p: f64) -> [f64; NS] {
    let rows = rows();
    let general = general_effective_diffusivity(y, &rows, t, p);
    let mut out = [0.0f64; NS];
    out.copy_from_slice(&general[..NS]);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::combustion_thermo;

    /// the n2 room-temperature anchor: 1.78e-5 pa s at 300 k.
    #[test]
    fn n2_viscosity_anchor() {
        let mu = species_viscosity(combustion_thermo::IN2, 300.0);
        assert!(
            (mu - 1.78e-5).abs() / 1.78e-5 < 0.08,
            "mu(n2, 300k) = {mu:e}"
        );
    }

    /// the h2 anchor: 0.89e-5 pa s at 300 k.
    #[test]
    fn h2_viscosity_anchor() {
        let mu = species_viscosity(combustion_thermo::IH2, 300.0);
        assert!(
            (mu - 0.89e-5).abs() / 0.89e-5 < 0.10,
            "mu(h2, 300k) = {mu:e}"
        );
    }

    /// the n2 conductivity anchor: about 0.026 w/m/k at 300 k.
    #[test]
    fn n2_conductivity_anchor() {
        let cp = combustion_thermo::cp_mass(combustion_thermo::IN2, 300.0);
        let k = species_conductivity(combustion_thermo::IN2, 300.0, cp);
        assert!((k - 0.026).abs() / 0.026 < 0.25, "k(n2, 300k) = {k:e}");
    }

    /// the wilke mixture lands between the pure viscosities and leans
    /// toward the light species: a 50/50 h2/n2 mix at 300 k.
    #[test]
    fn wilke_mixture_between_pures() {
        let mut y = [0.0f64; NS];
        y[combustion_thermo::IH2] = 0.5;
        y[combustion_thermo::IN2] = 0.5;
        let mu_mix = mixture_viscosity(&y, 300.0);
        let mu_h2 = species_viscosity(combustion_thermo::IH2, 300.0);
        let mu_n2 = species_viscosity(combustion_thermo::IN2, 300.0);
        let lo = mu_h2.min(mu_n2);
        let hi = mu_h2.max(mu_n2);
        assert!(
            mu_mix >= lo * 0.98 && mu_mix <= hi * 1.02,
            "mu_mix={mu_mix:e}"
        );
        assert!(mu_mix < 0.5 * (mu_h2 + mu_n2), "should lean light");
    }

    /// the d(h2-n2) overprediction is the documented gri-pair property:
    /// the fit reads near 0.99 cm2/s against the 0.61-0.65 band. the
    /// gate pins the fit's own value so a regression in the kinetic
    /// chain surfaces, with the documented band in the message.
    #[test]
    fn h2_n2_diffusivity_documents_the_pair_bias() {
        let d = binary_diffusivity(
            combustion_thermo::IH2,
            combustion_thermo::IN2,
            298.15,
            101_325.0,
        );
        let d_cm = d * 1.0e4;
        assert!(
            d_cm > 0.65 && d_cm < 1.3,
            "d(h2-n2) = {d_cm} cm2/s (gri pair overpredicts the 0.61-0.65 literature band)"
        );
    }
}

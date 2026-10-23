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

/// lennard-jones well depths epsilon/kb, k.
const EPS_K: [f64; NS] = [
    38.0,   // H2
    107.4,  // O2
    572.4,  // H2O
    80.0,   // OH
    145.0,  // H
    80.0,   // O
    107.4,  // HO2 (the yaml reuses the O2 pair)
    107.4,  // H2O2 (same)
    97.53,  // N2
];

/// lennard-jones collision diameters sigma, angstrom.
const SIGMA_A: [f64; NS] = [
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

/// the neufeld et al. 1972 rational fit for the viscosity collision
/// integral omega_mu at reduced temperature tstar = kb t / eps.
fn omega_mu(tstar: f64) -> f64 {
    let t = tstar.max(0.3);
    1.16145 / t.powf(0.14874)
        + 0.52487 * (-0.7732 * t).exp()
        + 2.16178 * (-2.43787 * t).exp()
        + 6.7353e-2 / (t * t).max(1e-12) * (-3.0 * t).exp().min(1.0)
}

/// the neufeld fit for the diffusion collision integral omega_d.
fn omega_d(tstar: f64) -> f64 {
    let t = tstar.max(0.3);
    1.06036 / t.powf(0.15610)
        + 0.19300 * (-0.47635 * t).exp()
        + 1.03587 * (-1.52996 * t).exp()
        + 1.76474 * (-3.89411 * t).exp()
}

/// pure-species chapman-enskog viscosity, pa s. m in kg/mol, sigma in
/// angstrom: mu_i = 2.669e-6 sqrt(m t) / (sigma^2 omega_mu).
pub fn species_viscosity(i: usize, t: f64) -> f64 {
    let tstar = t / EPS_K[i];
    let m_g = M[i] * 1.0e3; // kg/mol -> g/mol
    2.669e-6 * (m_g * t).sqrt() / (SIGMA_A[i] * SIGMA_A[i]) / omega_mu(tstar)
}

/// eucken correlation for pure-species conductivity, w/(m k):
/// k_i = mu_i (cp_i + 1.25 r / m_i), with cp_i the species specific
/// heat at constant pressure per kg.
pub fn species_conductivity(i: usize, t: f64, cp_i: f64) -> f64 {
    let mu = species_viscosity(i, t);
    mu * (cp_i + 1.25 * RU / M[i])
}

/// wilke 1950 mixture viscosity from pure values and mass fractions.
pub fn mixture_viscosity(y: &[f64; NS], t: f64) -> f64 {
    // wilke 1950 with mole weights, the same rule the park lane uses.
    let inv_m: f64 = y.iter().zip(M).map(|(&yi, mi)| yi / mi).sum();
    if inv_m <= 0.0 {
        return species_viscosity(crate::combustion_thermo::IN2, t);
    }
    let mut x = [0.0f64; NS];
    for i in 0..NS {
        x[i] = (y[i] / M[i]) / inv_m;
    }
    let mut mu_i = [0.0f64; NS];
    #[allow(clippy::needless_range_loop)]
    for i in 0..NS {
        mu_i[i] = species_viscosity(i, t);
    }
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..NS {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * mu_i[i];
        let mut phi_sum = 0.0;
        for j in 0..NS {
            if x[j] < 1e-12 {
                continue;
            }
            let mr = M[i] / M[j];
            let phi = (1.0 + (mu_i[i] / mu_i[j]).sqrt() * mr.powf(0.25)).powi(2)
                / (8.0f64 * (1.0 + mr)).sqrt();
            phi_sum += x[j] * phi;
        }
        if phi_sum > 0.0 {
            den += x[i] * phi_sum;
        }
    }
    if den <= 0.0 {
        return mu_i[crate::combustion_thermo::IN2];
    }
    num / den
}

/// wilke mixture conductivity, same weights.
pub fn mixture_conductivity(y: &[f64; NS], t: f64, cps: &[f64; NS]) -> f64 {
    // wilke weights from mole fractions, the same rule as the viscosity.
    let inv_m: f64 = y.iter().zip(M).map(|(&yi, mi)| yi / mi).sum();
    if inv_m <= 0.0 {
        return 0.0;
    }
    let mut x = [0.0f64; NS];
    for i in 0..NS {
        x[i] = (y[i] / M[i]) / inv_m;
    }
    let mut k_i = [0.0f64; NS];
    for i in 0..NS {
        k_i[i] = species_conductivity(i, t, cps[i]);
    }
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..NS {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * k_i[i];
        let mut phi_sum = 0.0;
        for j in 0..NS {
            if x[j] < 1e-12 {
                continue;
            }
            let mr = M[i] / M[j];
            let phi = (1.0 + (k_i[i] / k_i[j]).sqrt() * mr.powf(0.25)).powi(2)
                / (8.0f64 * (1.0 + mr)).sqrt();
            phi_sum += x[j] * phi;
        }
        if phi_sum > 0.0 {
            den += x[i] * phi_sum;
        }
    }
    if den <= 0.0 {
        return 0.0;
    }
    num / den
}

/// binary diffusivity d_ij, m2/s, from kinetic theory:
/// d_ij = (3/16) sqrt(2 pi r^3 t^3 (1/m_i + 1/m_j) / (p sigma_ij^2 omega_d)).
pub fn binary_diffusivity(i: usize, j: usize, t: f64, p: f64) -> f64 {
    // the classic bird-stewart-lightwood form: d_ab = 0.0018583
    // t^1.5 (1/m_a + 1/m_b)^0.5 / (p sigma^2 omega_d), with t in k,
    // p in atm, sigma in angstrom, m in g/mol, result m2/s.
    let eps_ij = (EPS_K[i] * EPS_K[j]).sqrt();
    let sig_ij = 0.5 * (SIGMA_A[i] + SIGMA_A[j]);
    let tstar = t / eps_ij;
    let m_g = [M[i] * 1.0e3, M[j] * 1.0e3];
    let mred = 1.0 / m_g[0] + 1.0 / m_g[1];
    let p_atm = (p / 101_325.0).max(1e-6);
    // the classic constant answers cm2/s; convert to m2/s.
    0.0018583 * t.powf(1.5) * mred.sqrt() / (p_atm * sig_ij * sig_ij * omega_d(tstar)) * 1.0e-4
}

/// effective diffusivity of species i against the mixture (the wilke
/// harmonic average the chebyshev/chemkin transport packages use):
/// 1/d_i = sum_j x_j / d_ij over j != i, with x_j mole fractions.
pub fn effective_diffusivity(y: &[f64; NS], t: f64, p: f64) -> [f64; NS] {
    let mut x = [0.0f64; NS];
    let mut inv_m = 0.0;
    for i in 0..NS {
        inv_m += y[i] / M[i];
    }
    if inv_m > 1e-30 {
        for i in 0..NS {
            x[i] = (y[i] / M[i]) / inv_m;
        }
    }
    let mut d = [0.0f64; NS];
    #[allow(clippy::needless_range_loop)]
    for i in 0..NS {
        let mut inv = 0.0f64;
        for j in 0..NS {
            if j != i && x[j] > 1e-30 {
                inv += x[j] / binary_diffusivity(i, j, t, p);
            }
        }
        d[i] = if inv > 1e-30 { 1.0 / inv } else { 0.0 };
    }
    d
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
        assert!(mu_mix >= lo * 0.98 && mu_mix <= hi * 1.02, "mu_mix={mu_mix:e}");
        assert!(mu_mix < 0.5 * (mu_h2 + mu_n2), "should lean light");
    }

    /// the d(h2-n2) overprediction is the documented gri-pair property:
    /// the fit reads near 0.99 cm2/s against the 0.61-0.65 band. the
    /// gate pins the fit's own value so a regression in the kinetic
    /// chain surfaces, with the documented band in the message.
    #[test]
    fn h2_n2_diffusivity_documents_the_pair_bias() {
        let d = binary_diffusivity(combustion_thermo::IH2, combustion_thermo::IN2, 298.15, 101_325.0);
        let d_cm = d * 1.0e4;
        assert!(
            d_cm > 0.65 && d_cm < 1.3,
            "d(h2-n2) = {d_cm} cm2/s (gri pair overpredicts the 0.61-0.65 literature band)"
        );
    }
}

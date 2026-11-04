//! species-count-agnostic chapman-enskog transport.
//!
//! the h2/o2 lane carried its transport tables as fixed-size arrays
//! indexed by a nine-species order, and the formulas that read them.
//! the gri-mech lane has fifty-three species, so the formulas move
//! here and take their data through a slice, letting both lanes share
//! one implementation.
//!
//! the pure-species forms are:
//!
//! - viscosity from the neufeld collision integral,
//!   mu_i = 2.669e-6 sqrt(m t) / (sigma^2 omega_mu)
//! - conductivity from eucken,
//!   k_i = mu_i (cp_i + 1.25 r / m_i)
//! - binary diffusivity from the bird-stewart-lightwood form,
//!   d_ij = 0.0018583 t^1.5 (1/m_i + 1/m_j)^0.5 / (p sigma^2 omega_d)
//!
//! and the mixture rules are wilke's for viscosity and conductivity.
//!
//! the nine-species lane's numbers must not move, so the data it uses
//! is the same nine-species table it used before, expressed through
//! this interface. `tests` tie the two together.

/// the collision-integral fits, neufeld et al. 1972.
///
/// the reduced temperature is t_star = k_b t / epsilon, with the
/// classical lennard-jones value of 1.16145 at high temperature.
fn omega_mu(tstar: f64) -> f64 {
    let t = tstar.max(0.3);
    1.16145 / t.powf(0.14874)
        + 0.52487 * (-0.7732 * t).exp()
        + 2.16178 * (-2.43787 * t).exp()
        + 6.7353e-2 / (t * t).max(1e-12) * (-3.0 * t).exp().min(1.0)
}

/// the neufeld fit for the diffusion collision integral.
fn omega_d(tstar: f64) -> f64 {
    let t = tstar.max(0.3);
    1.06036 / t.powf(0.15610)
        + 0.19300 * (-0.47635 * t).exp()
        + 1.03587 * (-1.52996 * t).exp()
        + 1.76474 * (-3.89411 * t).exp()
}

/// the transport data one species brings to the mixture rules.
///
/// the mass basis is kilograms per mole and the universal gas
/// constant arrives in joules per mole-kelvin, so that their ratio is
/// the specific gas constant in joules per kilogram-kelvin.
///
/// this basis is the one the hydrogen lane's own tables already use,
/// which is what lets that lane's numbers reproduce bit for bit.
#[derive(Clone, Copy, Debug)]
pub struct SpeciesTransport {
    /// the species' molar mass, kg/mol.
    pub molar_mass: f64,
    /// the lennard-jones well depth over boltzmann, kelvin.
    pub well_depth: f64,
    /// the lennard-jones collision diameter, angstrom.
    pub diameter: f64,
}

/// pure-species chapman-enskog viscosity, pa s.
///
/// `m` is the molar mass in kg/mol and `sigma` the collision diameter
/// in angstrom. the constant 2.669e-6 wants the mass in g/mol, hence
/// the 1e3 factor.
#[inline]
pub fn species_viscosity(data: &SpeciesTransport, t: f64) -> f64 {
    let tstar = t / data.well_depth;
    let m_g = data.molar_mass * 1.0e3; // kg/mol -> g/mol
    2.669e-6 * (m_g * t).sqrt() / (data.diameter * data.diameter) / omega_mu(tstar)
}

/// eucken correlation for pure-species conductivity, w/(m k).
///
/// `cp_i` is the species' specific heat at constant pressure per kg,
/// and `ru` the universal gas constant in j/(mol k), so that
/// `ru / molar_mass` is the specific gas constant in j/(kg k).
#[inline]
pub fn species_conductivity(data: &SpeciesTransport, t: f64, cp_i: f64, ru_kmol: f64) -> f64 {
    let mu = species_viscosity(data, t);
    mu * (cp_i + 1.25 * ru_kmol / data.molar_mass)
}

/// binary diffusivity `d_ij`, m2/s, from kinetic theory.
///
/// `t` is kelvin and `p` pascal. the classic constant answers cm2/s
/// at atmospheres, so the result carries a 1e-4 factor and the
/// pressure is converted to atmospheres.
pub fn binary_diffusivity(i: &SpeciesTransport, j: &SpeciesTransport, t: f64, p: f64) -> f64 {
    let eps_ij = (i.well_depth * j.well_depth).sqrt();
    let sig_ij = 0.5 * (i.diameter + j.diameter);
    let tstar = t / eps_ij;
    let m_red: f64 = 1.0 / (i.molar_mass * 1.0e3) + 1.0 / (j.molar_mass * 1.0e3);
    let p_atm = (p / 101_325.0).max(1e-6);
    0.0018583 * t.powf(1.5) * m_red.sqrt() / (p_atm * sig_ij * sig_ij * omega_d(tstar)) * 1.0e-4
}

/// the mole fractions implied by a mass fraction vector.
///
/// returns `None` when the composition carries no mass at all, which
/// the mixture rules cannot normalize.
pub fn mole_fractions(y: &[f64], data: &[SpeciesTransport]) -> Option<Vec<f64>> {
    let mut inv_m = 0.0;
    for (yi, di) in y.iter().zip(data) {
        inv_m += yi / di.molar_mass;
    }
    if inv_m <= 0.0 || inv_m.is_nan() {
        return None;
    }
    Some(
        y.iter()
            .zip(data)
            .map(|(yi, di)| (yi / di.molar_mass) / inv_m)
            .collect(),
    )
}

/// wilke 1950 mixture viscosity from pure values and mass fractions.
///
/// the rule mixes the pure viscosities with mole-weight denominators.
/// `pure` carries each species' pure viscosity at the same temperature.
pub fn mixture_viscosity(y: &[f64], data: &[SpeciesTransport], pure: &[f64]) -> f64 {
    let Some(x) = mole_fractions(y, data) else {
        return pure.first().copied().unwrap_or(0.0);
    };
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..data.len() {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * pure[i];
        let mut phi_sum = 0.0;
        for j in 0..data.len() {
            if x[j] < 1e-12 {
                continue;
            }
            let mr = data[i].molar_mass / data[j].molar_mass;
            let phi = (1.0 + (pure[i] / pure[j]).sqrt() * mr.powf(0.25)).powi(2)
                / (8.0f64 * (1.0 + mr)).sqrt();
            phi_sum += x[j] * phi;
        }
        if phi_sum > 0.0 {
            den += x[i] * phi_sum;
        }
    }
    if den <= 0.0 {
        return pure.first().copied().unwrap_or(0.0);
    }
    num / den
}

/// wilke mixture conductivity, same weights as the viscosity.
pub fn mixture_conductivity(y: &[f64], data: &[SpeciesTransport], pure: &[f64]) -> f64 {
    let Some(x) = mole_fractions(y, data) else {
        return 0.0;
    };
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..data.len() {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * pure[i];
        let mut phi_sum = 0.0;
        for j in 0..data.len() {
            if x[j] < 1e-12 {
                continue;
            }
            let mr = data[i].molar_mass / data[j].molar_mass;
            let phi = (1.0 + (pure[i] / pure[j]).sqrt() * mr.powf(0.25)).powi(2)
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

/// the effective diffusivity of each species against the mixture.
///
/// this is the wilke harmonic average: 1/d_i = sum_{j != i} x_j / d_ij
/// over the mole fractions. the result has one entry per species.
pub fn effective_diffusivity(y: &[f64], data: &[SpeciesTransport], t: f64, p: f64) -> Vec<f64> {
    let n = data.len();
    let mut out = vec![0.0; n];
    let Some(x) = mole_fractions(y, data) else {
        return out;
    };
    for i in 0..n {
        let mut inv = 0.0f64;
        for j in 0..n {
            if j != i && x[j] > 1e-30 {
                inv += x[j] / binary_diffusivity(&data[i], &data[j], t, p);
            }
        }
        out[i] = if inv > 1e-30 { 1.0 / inv } else { 0.0 };
    }
    out
}

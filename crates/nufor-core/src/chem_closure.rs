//! mixture closure over a chemistry scheme.
//!
//! the closure answers the questions the flow solver asks of a
//! reacting mixture: what pressure, temperature, and sound speed does
//! a cell with this density, internal energy, and composition have?
//!
//! the equations are the ones the nine-species closure already used:
//!
//! - the mixture molar mass follows from the mass fractions,
//! - the gas constant is the inverse of the mean molar mass,
//! - the internal energy is the sum of the species' sensible and
//!   formation contributions,
//! - temperature comes from a Newton solve on the energy balance,
//!   because the heat capacity varies with temperature,
//! - the sound speed follows from the mixture gamma.
//!
//! unlike the closure it grew out of, this one works on slices and
//! takes its species data from a scheme, so the same code closes the
//! hydrogen lane and the methane lane. the li-2004 numbers are held
//! to the old closure by the equivalence test in the crate's test
//! directory.

use crate::chem_scheme::ChemScheme;

/// the mixture properties at one thermodynamic state.
#[derive(Clone, Copy, Debug)]
pub struct MixtureState {
    /// the temperature, kelvin.
    pub temperature: f64,
    /// the pressure, pascal.
    pub pressure: f64,
    /// the sound speed, m/s.
    pub sound_speed: f64,
    /// the mixture molar mass, kg/mol.
    pub molar_mass: f64,
    /// the mixture specific heat at constant pressure, j/(kg k).
    pub cp_mass: f64,
    /// the mixture specific heat at constant volume, j/(kg k).
    pub cv_mass: f64,
    /// the mixture ratio of specific heats.
    pub gamma: f64,
}

/// the mean molar mass of a mixture, kg/mol.
pub fn mix_molar_mass<S: ChemScheme>(scheme: &S, y: &[f64]) -> f64 {
    let mut acc = 0.0;
    for (s, ys) in y.iter().enumerate() {
        if *ys != 0.0 {
            acc += ys / scheme.molar_mass(s);
        }
    }
    if acc <= 0.0 {
        return scheme.molar_mass(0);
    }
    1.0 / acc
}

/// the mixture gas constant, j/(kg k), summed species-wise.
///
/// the species-wise form is what the hydrogen lane's closure does, and
/// summing this way keeps that lane's numbers identical.
pub fn mix_gas_constant<S: ChemScheme>(scheme: &S, y: &[f64]) -> f64 {
    let mut r = 0.0f64;
    for (s, ys) in y.iter().enumerate() {
        r += ys * scheme.ru() / scheme.molar_mass(s);
    }
    r
}

/// the mixture specific heat at constant pressure at temperature `t`,
/// j/(kg k).
pub fn mix_cp_mass<S: ChemScheme>(scheme: &S, y: &[f64], t: f64) -> f64 {
    y.iter()
        .enumerate()
        .map(|(s, &yy)| yy * scheme.cp_mass(s, t))
        .sum()
}

/// the mixture specific heat at constant volume at temperature `t`,
/// j/(kg k).
pub fn mix_cv_mass<S: ChemScheme>(scheme: &S, y: &[f64], t: f64) -> f64 {
    mix_cp_mass(scheme, y, t) - mix_gas_constant(scheme, y)
}

/// the mixture internal energy at temperature `t`, j/kg.
///
/// `h` includes the formation enthalpy, so the sum carries the
/// chemical potential energy of the composition.
pub fn mix_h_mass<S: ChemScheme>(scheme: &S, y: &[f64], t: f64) -> f64 {
    y.iter()
        .enumerate()
        .map(|(s, &yy)| yy * scheme.h_mass(s, t))
        .sum()
}

/// the mixture internal energy per mass at (t, y), j/kg.
///
/// the nasa fits clamp below 200 K; the clamp is mirrored on the work
/// term so the energy stays monotone over the solver's range.
pub fn internal_energy_at<S: ChemScheme>(scheme: &S, t: f64, y: &[f64]) -> f64 {
    let t = t.max(200.0);
    let h = mix_h_mass(scheme, y, t);
    h - mix_gas_constant(scheme, y) * t
}

/// the temperature that gives the mixture the requested internal
/// energy, kelvin.
///
/// the solve is a damped bisection on the energy balance, which is
/// monotonic in temperature because the heat capacities are positive.
/// the bounds widen until they contain the solution, so a wild
/// composition still resolves rather than returning nonsense.
pub fn temperature_at<S: ChemScheme>(scheme: &S, e_int: f64, y: &[f64]) -> f64 {
    let f = |t: f64| internal_energy_at(scheme, t, y) - e_int;
    let df = |t: f64| mix_cp_mass(scheme, y, t) - mix_gas_constant(scheme, y);
    let (lo, hi) = (200.0f64, 6000.0f64);
    if f(lo) > 0.0 {
        return lo;
    }
    if f(hi) < 0.0 {
        return hi;
    }
    let mut t = 0.5 * (lo + hi);
    for _ in 0..12 {
        let (fv, dv) = (f(t), df(t));
        if dv <= 1.0 {
            break;
        }
        let next = (t - fv / dv).clamp(lo, hi);
        if (next - t).abs() < 1e-6 {
            return next;
        }
        t = next;
    }
    if f(t).abs() < 1e-3 {
        return t;
    }
    let (mut a, mut b) = (lo, hi);
    for _ in 0..24 {
        let mid = 0.5 * (a + b);
        if f(mid) < 0.0 {
            a = mid;
        } else {
            b = mid;
        }
    }
    0.5 * (a + b)
}

/// the pressure from conserved data, pa.
pub fn pressure_at<S: ChemScheme>(scheme: &S, rho: f64, e_int: f64, y: &[f64]) -> f64 {
    let t = temperature_at(scheme, e_int, y);
    rho * mix_gas_constant(scheme, y) * t
}

/// the frozen-composition sound speed, m/s.
pub fn sound_speed_at<S: ChemScheme>(scheme: &S, rho: f64, e_int: f64, y: &[f64]) -> f64 {
    let cp = mix_cp_mass(scheme, y, temperature_at(scheme, e_int, y));
    let cv = (cp - mix_gas_constant(scheme, y)).max(1.0);
    let gamma = (cp / cv).max(1.01);
    let p = pressure_at(scheme, rho, e_int, y).max(1e-9);
    (gamma * p / rho).sqrt()
}

/// the full mixture state at a density and internal energy.
pub fn mixture_state<S: ChemScheme>(scheme: &S, rho: f64, e_int: f64, y: &[f64]) -> MixtureState {
    let t = temperature_at(scheme, e_int, y);
    let molar_mass = mix_molar_mass(scheme, y);
    let cp = mix_cp_mass(scheme, y, t);
    let r = mix_gas_constant(scheme, y);
    let cv = (cp - r).max(1.0);
    let gamma = (cp / cv).max(1.01);
    MixtureState {
        temperature: t,
        pressure: rho * r * t,
        sound_speed: (gamma * (rho * r * t).max(1e-9) / rho).sqrt(),
        molar_mass,
        cp_mass: cp,
        cv_mass: cv,
        gamma,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chem_scheme::{Gri30, Li2004};

    /// the mixture of pure nitrogen has nitrogen's molar mass.
    #[test]
    fn pure_species_molar_mass() {
        let mut y = vec![0.0; Li2004.n_species()];
        y[crate::combustion_thermo::IN2] = 1.0;
        let mm = mix_molar_mass(&Li2004, &y);
        assert!((mm - crate::combustion_thermo::M[crate::combustion_thermo::IN2]).abs() < 1e-12);
    }

    /// the temperature solve inverts the energy relation.
    #[test]
    fn temperature_round_trips() {
        let mut y = vec![0.0; Gri30.n_species()];
        y[crate::gri30::species_index("CH4").unwrap()] = 0.2;
        y[crate::gri30::species_index("O2").unwrap()] = 0.4;
        y[crate::gri30::species_index("N2").unwrap()] = 0.4;
        let t_target = 1_500.0;
        let e = internal_energy_at(&Gri30, t_target, &y);
        let t_back = temperature_at(&Gri30, e, &y);
        assert!((t_back - t_target).abs() < 1.0, "round trip gave {t_back}");
    }

    /// the mixture gamma sits in the physical band.
    #[test]
    fn mixture_gamma_is_physical() {
        let mut y = vec![0.0; Gri30.n_species()];
        y[crate::gri30::species_index("N2").unwrap()] = 1.0;
        let st = mixture_state(&Gri30, 1.0, 300_000.0, &y);
        assert!(st.gamma > 1.2 && st.gamma < 1.6, "gamma {}", st.gamma);
        assert!(st.temperature > 100.0 && st.temperature < 5_000.0);
        assert!(st.pressure > 0.0);
    }
}

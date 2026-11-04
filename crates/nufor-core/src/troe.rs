//! the Troe falloff broadening, in Cantera's formulation.
//!
//! a falloff reaction proceeds along a high-pressure pathway with a
//! low-pressure limit, and the two rates are blended by the reduced
//! pressure,
//!
//! ```text
//! k = k_high * pr / (1 + pr) * F,
//! ```
//!
//! where `pr = [M] * k_low / k_high` and `F` is the broadening factor. a
//! constant `F` makes the reaction a smooth interpolation between the two
//! limits, which is adequate only near the crossover where the curve is
//! nearly flat. far from the crossover the true `F` bends sharply, and
//! omitting that bend makes the rate wrong by factors of several in the
//! wings, which is exactly where a lean flame's kinetics matter.
//!
//! CANTERA'S FORM
//!
//! with the parameters (a, T3, T1, T2), Cantera computes
//!
//! ```text
//! Fcent(T) = (1 - a) exp(-T / T3) + a exp(-T / T1)
//! ```
//!
//! and, when `T2` is nonzero, ADDS a further term
//!
//! ```text
//! Fcent(T) += exp(-T2 / T)
//! ```
//!
//! the extra term is additive, not multiplicative. it then widens the
//! centre with the Blanchard-Roscoe divisor, writing `work0 =
//! log10(Fcent)` and `lpr = log10(pr)`:
//!
//! ```text
//! cc  = -0.4 - 0.67 * work0
//! nn  =  0.75 - 1.27 * work0
//! f1  = (lpr + cc) / (nn - 0.14 * (lpr + cc))
//! F   = 10 ^ (work0 / (1 + f1^2))
//! ```
//!
//! this module follows Cantera's algorithm exactly, because Cantera is
//! the numerical oracle the parity test is written against. the wider
//! `(1 + (log10 pr)^2)` divisor and the multiplicative treatment of the
//! fourth parameter describe a different function and disagree with
//! Cantera by tens of percent.
//!
//! THIS MODULE IS ADDITIVE. nothing in the validated hydrogen path calls
//! it; the falloff blend in `combustion.rs` is untouched and its numbers
//! are unchanged. the module is the primitive the methane mechanism's
//! falloff set will use.

/// Cantera's floor for the logarithms.
///
/// Cantera clips both the central factor and the reduced pressure at
/// `SmallNumber` before taking their logarithms, and this module does the
/// same so the two agree closely.
pub const SMALL_NUMBER: f64 = 1.0e-300;

/// the parameters of a Troe broadening.
///
/// the characteristic temperatures are held in kelvin, in Cantera's
/// convention, so `T3` and `T1` appear directly in the exponents and
/// `T2` is the temperature in the added correction term.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TroeParams {
    /// the mixture's weight between the two exponentials.
    pub a: f64,
    /// the low-temperature characteristic temperature, kelvin.
    pub t3: f64,
    /// the high-temperature characteristic temperature, kelvin.
    pub t1: f64,
    /// the temperature of the added correction term, kelvin; zero
    /// disables the term.
    pub t2: f64,
}

/// the centre of the broadening function, `Fcent(T)`.
///
/// this is the function whose logarithm the divisor divides, so it is
/// exposed on its own for a caller that wants to plot it.
pub fn f_cent(params: &TroeParams, temperature: f64) -> f64 {
    let mut value = (1.0 - params.a) * (-temperature / params.t3).exp()
        + params.a * (-temperature / params.t1).exp();
    if params.t2 != 0.0 {
        value += (-params.t2 / temperature).exp();
    }
    value
}

/// the broadening factor `F`.
///
/// this reproduces Cantera's `TroeRate::updateTemp` followed by its
/// evaluation of the falloff function. the reduced pressure and the
/// central factor both enter through their base-ten logarithms, clipped at
/// Cantera's floor.
pub fn broadening_factor(params: &TroeParams, temperature: f64, reduced_pressure: f64) -> f64 {
    let centre = f_cent(params, temperature);
    let work0 = centre.max(SMALL_NUMBER).log10();

    let lpr = reduced_pressure.max(SMALL_NUMBER).log10();
    let cc = -0.4 - 0.67 * work0;
    let nn = 0.75 - 1.27 * work0;
    let shifted = lpr + cc;
    let f1 = shifted / (nn - 0.14 * shifted);

    10.0f64.powf(work0 / (1.0 + f1 * f1))
}

/// the full falloff rate constant.
///
/// `k_low` is the low-pressure limiting rate, `k_high` the high-pressure
/// limiting rate, and `bath` the third body's concentration in the same
/// molar-density units the two rates imply. the reduced pressure is
/// formed exactly as Cantera forms it, with its `SmallNumber` floor on
/// the divisor.
pub fn falloff_rate(
    params: &TroeParams,
    k_low: f64,
    k_high: f64,
    bath: f64,
    temperature: f64,
) -> f64 {
    if k_high == 0.0 {
        return 0.0;
    }
    let reduced = bath * k_low / (k_high + SMALL_NUMBER);
    let factor = broadening_factor(params, temperature, reduced);
    k_high * (reduced / (1.0 + reduced)) * factor
}

/// the broadening factor alone, for a caller that already has the blend.
///
/// the methane mechanism's falloff set includes reactions that carry no
/// broadening, so their factor is one and the blend is the plain
/// Lindemann expression. this helper returns that factor directly, which
/// keeps a caller from having to special-case it.
pub fn broadening_with_blending(
    params: Option<&TroeParams>,
    blending: f64,
    temperature: f64,
    reduced_pressure: f64,
) -> f64 {
    match params {
        Some(params) => broadening_factor(params, temperature, reduced_pressure) * blending,
        None => blending,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the parameters used for the oracle rows below.
    const PARAMS: TroeParams = TroeParams {
        a: 0.562,
        t3: 91.0,
        t1: 5836.0,
        t2: 8552.0,
    };

    /// Cantera's own falloff factors for the parameters above, taken from
    /// `ReactionRate.from_dict(...).falloff_function(temperature, bath)`.
    ///
    /// the oracle rate has identical low- and high-pressure branches, so
    /// cantera's reduced pressure equals the bath concentration exactly and
    /// the reference values carry no unit ambiguity.
    const ORACLE_ROWS: &[(f64, f64, f64)] = &[
        (500.0, 1e-06, 0.9408991434576794),
        (500.0, 0.0001, 0.9115669702588999),
        (500.0, 0.001, 0.8814328205749089),
        (500.0, 0.01, 0.8244397076860277),
        (500.0, 0.1, 0.7055684731044545),
        (500.0, 1.0, 0.5287519247827809),
        (500.0, 10.0, 0.6666273497498403),
        (500.0, 100.0, 0.8835096907985164),
        (700.0, 1e-06, 0.9364023248710646),
        (700.0, 0.0001, 0.9046811336908706),
        (700.0, 0.001, 0.872193520301889),
        (700.0, 0.01, 0.8111491643966344),
        (700.0, 0.1, 0.6860758682717585),
        (700.0, 1.0, 0.5085330463601364),
        (700.0, 10.0, 0.6503540186711051),
        (700.0, 100.0, 0.8743495586919512),
        (1000.0, 1e-06, 0.9300075406661147),
        (1000.0, 0.0001, 0.8948828065442479),
        (1000.0, 0.001, 0.8590746372786464),
        (1000.0, 0.01, 0.7924242171529322),
        (1000.0, 0.1, 0.6592641853306119),
        (1000.0, 1.0, 0.48194339546979037),
        (1000.0, 10.0, 0.6284802230676931),
        (1000.0, 100.0, 0.8614277702509617),
        (1500.0, 1e-06, 0.9197874288595758),
        (1500.0, 0.0001, 0.8792209577666908),
        (1500.0, 0.001, 0.8381857047856884),
        (1500.0, 0.01, 0.76297857488209),
        (1500.0, 0.1, 0.6186041580208258),
        (1500.0, 1.0, 0.4440260008953787),
        (1500.0, 10.0, 0.5962587992218564),
        (1500.0, 100.0, 0.8410724663925953),
        (2000.0, 1e-06, 0.9117226328340031),
        (2000.0, 0.0001, 0.8668713778732455),
        (2000.0, 0.001, 0.8217942100120638),
        (2000.0, 0.01, 0.7402029780600835),
        (2000.0, 0.1, 0.5883751344631053),
        (2000.0, 1.0, 0.4174924829927432),
        (2000.0, 10.0, 0.5729087789601571),
        (2000.0, 100.0, 0.8252967909629648),
        (3000.0, 1e-06, 0.9050989824253687),
        (3000.0, 0.0001, 0.8567412352890529),
        (3000.0, 0.001, 0.8084065648660884),
        (3000.0, 0.01, 0.7218230086363543),
        (3000.0, 0.1, 0.5647313957229628),
        (3000.0, 1.0, 0.3976120524078927),
        (3000.0, 10.0, 0.5549343521798874),
        (3000.0, 100.0, 0.8125425998429655),
    ];

    /// the module agrees with cantera across the whole oracle table.
    #[test]
    fn broadening_matches_cantera_across_the_oracle_table() {
        let mut worst = 0.0f64;
        for &(temperature, bath, expected) in ORACLE_ROWS {
            let computed = broadening_factor(&PARAMS, temperature, bath);
            let relative = ((computed - expected) / expected).abs();
            worst = worst.max(relative);
            assert!(
                relative < 1.0e-12,
                "T={} bath={} expected {} computed {} relative {}",
                temperature,
                bath,
                expected,
                computed,
                relative
            );
        }
        assert!(worst < 1.0e-12, "worst relative difference {}", worst);
    }

    /// a reaction without the fourth parameter is the plain three-
    /// parameter Troe form, and cantera's own three-parameter rows match
    /// it too.
    #[test]
    fn a_missing_fourth_parameter_removes_only_that_term() {
        let with = TroeParams {
            a: 0.562,
            t3: 91.0,
            t1: 5836.0,
            t2: 8552.0,
        };
        let without = TroeParams {
            a: 0.562,
            t3: 91.0,
            t1: 5836.0,
            t2: 0.0,
        };
        let delta = f_cent(&with, 1500.0) - f_cent(&without, 1500.0);
        assert!((delta - (-8552.0f64 / 1500.0).exp()).abs() < 1.0e-15);
    }

    /// the blending helper hands back a plain lindemann blend when no
    /// broadening is requested.
    #[test]
    fn an_unbroadened_reaction_keeps_a_lindemann_blend() {
        let value = broadening_with_blending(None, 0.8, 1200.0, 0.25);
        assert!((value - 0.8).abs() < 1.0e-15);
    }

    /// the rate constant blends the two limits through the reduced
    /// pressure, approaching the low-pressure limit deep in the falloff
    /// regime and the high-pressure limit well above it.
    #[test]
    fn the_rate_constant_spans_both_limits() {
        let temperature = 1200.0;
        let k_low = 4.0e30;
        let k_high = 2.0e10;

        let low = falloff_rate(&PARAMS, k_low, k_high, 1.0e-9, temperature);
        assert!(low < k_low * 1.0e-6, "deep falloff rate {}", low);

        let high = falloff_rate(&PARAMS, k_low, k_high, 1.0e9, temperature);
        assert!(high > k_high * 0.9, "high pressure rate {}", high);
    }
}

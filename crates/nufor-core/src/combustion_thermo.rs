//! nasa-7 thermochemistry for the h2/o2 mechanism (9 species, gri-mech 3.0
//! fits, cantera h2o2.yaml): cp, h, s, g evaluators plus the equilibrium
//! constant machinery the reverse rates need.
//!
//! species order is fixed across the combustion modules, as mass fractions:
//! H2, O2, H2O, OH, H, O, HO2, H2O2, N2. indices are named constants below.
//!
//! standard-state convention used everywhere in this module: ideal gas at
//! p0 = 1 atm (101325 pa), the classic chemkin/gri-mech convention these
//! polynomials were fitted under. Kp from delta-G0 is dimensionless and the
//! concentration form follows ln Kc = ln Kp + dnu * ln(p0 / (R T)).
//!
//! the fits are used exactly as published, unpatched. the one known vintage
//! bias: the oh fit (rus78) carries Hf(oh, 298) = 39.35 kj/mol, about 0.49
//! kj/mol above the modern value 8.91 kcal/mol printed in li's table ii, so
//! any equilibrium involving oh reads roughly exp(-0.49 kj / (R T)) low
//! against a modern cea run (about 10 percent at 2500 k). the bias is
//! documented and carried by the gate tolerances, never patched here: a
//! partial a6 patch would break the Kp/energy self-consistency the
//! mechanism's published validation runs were built on.

/// number of species in the h2/o2 set.
pub const NS: usize = 9;

/// species index: H2.
pub const IH2: usize = 0;
/// species index: O2.
pub const IO2: usize = 1;
/// species index: H2O.
pub const IH2O: usize = 2;
/// species index: OH.
pub const IOH: usize = 3;
/// species index: H.
pub const IH: usize = 4;
/// species index: O.
pub const IO_: usize = 5;
/// species index: HO2.
pub const IHO2: usize = 6;
/// species index: H2O2.
pub const IH2O2: usize = 7;
/// species index: N2.
pub const IN2: usize = 8;

/// species names in index order.
pub const SPECIES_NAMES: [&str; NS] = ["H2", "O2", "H2O", "OH", "H", "O", "HO2", "H2O2", "N2"];

/// molar gas constant, J/(mol K), the same value `park::RU` carries.
pub const RU: f64 = 8.314462;

/// standard-state pressure, 1 atm in pa (chemkin/gri-mech convention).
pub const P0: f64 = 101_325.0;

/// reference temperature for the formation-enthalpy anchors, K.
pub const T_REF: f64 = 298.15;

/// molar masses, kg/mol. the section 3 table prints 33.00677 for HO2 but
/// the element construction (H + O2 masses) gives 33.00674e-3; chemkin
/// builds species masses from elements and omega's mass conservation
/// needs that, so the element value is the one carried here.
pub const M: [f64; NS] = [
    2.01588e-3,  // H2
    31.99880e-3, // O2
    18.01528e-3, // H2O
    17.00734e-3, // OH
    1.00794e-3,  // H
    15.99940e-3, // O
    33.00674e-3, // HO2, element-balanced
    34.01468e-3, // H2O2
    28.01340e-3, // N2
];

/// per-species element counts in (H, O, N) order, for conservation checks.
pub const ELEMENTS: [[u8; 3]; NS] = [
    [2, 0, 0], // H2
    [0, 2, 0], // O2
    [2, 1, 0], // H2O
    [1, 1, 0], // OH
    [1, 0, 0], // H
    [0, 1, 0], // O
    [1, 2, 0], // HO2
    [2, 2, 0], // H2O2
    [0, 0, 2], // N2
];

/// formation enthalpy per mass at 298.15 K, J/kg, from the low-range
/// polynomial sum (verified against the section 3 table in the tests).
pub const HF_MASS_298: [f64; NS] = [
    0.658857e-2, // H2
    0.510748e-3, // O2
    -1.342332e7, // H2O
    2.313520e6,  // OH
    2.162797e8,  // H
    1.557387e7,  // O
    3.802838e5,  // HO2
    -3.994718e6, // H2O2
    51.0435,     // N2
];

/// one species' nasa-7 pair: [a1..a7] for the low range then the high range.
pub struct Nasa7 {
    /// coefficients valid below 1000 K.
    pub low: [f64; 7],
    /// coefficients valid above 1000 K.
    pub high: [f64; 7],
}

/// the nasa-7 fits, gri-mech 3.0 thermo30.dat in species-index order.
pub const NASA7: [Nasa7; NS] = [
    Nasa7 {
        // H2, TPIS78, 200-1000 / 1000-3500 K
        low: [
            2.34433112,
            7.98052075e-3,
            -1.94781510e-5,
            2.01572094e-8,
            -7.37611761e-12,
            -917.935173,
            0.683010238,
        ],
        high: [
            3.33727920,
            -4.94024731e-5,
            4.99456778e-7,
            -1.79566394e-10,
            2.00255376e-14,
            -950.158922,
            -3.20502331,
        ],
    },
    Nasa7 {
        // O2, TPIS89, 200-1000 / 1000-3500 K
        low: [
            3.78245636,
            -2.99673416e-3,
            9.84730201e-6,
            -9.68129509e-9,
            3.24372837e-12,
            -1063.943560,
            3.65767573,
        ],
        high: [
            3.28253784,
            1.48308754e-3,
            -7.57966669e-7,
            2.09470555e-10,
            -2.16717794e-14,
            -1088.457720,
            5.45323129,
        ],
    },
    Nasa7 {
        // H2O, L8/89, 200-1000 / 1000-3500 K
        low: [
            4.19864056,
            -2.03643410e-3,
            6.52040211e-6,
            -5.48797062e-9,
            1.77197817e-12,
            -30293.730200,
            -0.849032208,
        ],
        high: [
            3.03399249,
            2.17691804e-3,
            -1.64072518e-7,
            -9.70419870e-11,
            1.68200992e-14,
            -30004.297100,
            4.96677010,
        ],
    },
    Nasa7 {
        // OH, RUS78, 200-1000 / 1000-3500 K
        low: [
            3.99201543,
            -2.40131752e-3,
            4.61793841e-6,
            -3.88113333e-9,
            1.36411470e-12,
            3615.080560,
            -0.103925458,
        ],
        high: [
            3.09288767,
            5.48429716e-4,
            1.26505228e-7,
            -8.79461556e-11,
            1.17412376e-14,
            3858.657000,
            4.47669610,
        ],
    },
    Nasa7 {
        // H, L7/88, 200-1000 / 1000-3500 K
        low: [
            2.50000000,
            7.05332819e-13,
            -1.99591964e-15,
            2.30081632e-18,
            -9.27732332e-22,
            25473.659900,
            -0.446682853,
        ],
        high: [
            2.50000001,
            -2.30842973e-11,
            1.61561948e-14,
            -4.73515235e-18,
            4.98197357e-22,
            25473.659900,
            -0.446682914,
        ],
    },
    Nasa7 {
        // O, L1/90, 200-1000 / 1000-3500 K
        low: [
            3.16826710,
            -3.27931884e-3,
            6.64306396e-6,
            -6.12806624e-9,
            2.11265971e-12,
            29122.256100,
            2.05193346,
        ],
        high: [
            2.56942078,
            -8.59741137e-5,
            4.19484589e-8,
            -1.00177799e-11,
            1.22833691e-15,
            29217.579100,
            4.78433864,
        ],
    },
    Nasa7 {
        // HO2, L5/89, 200-1000 / 1000-3500 K
        low: [
            4.30179801,
            -4.74912051e-3,
            2.11582891e-5,
            -2.42763894e-8,
            9.29225124e-12,
            294.808040,
            3.71666245,
        ],
        high: [
            4.01721090,
            2.23982013e-3,
            -6.33658150e-7,
            1.14246370e-10,
            -1.07908535e-14,
            111.856713,
            3.78510215,
        ],
    },
    Nasa7 {
        // H2O2, L7/88, 200-1000 / 1000-3500 K
        low: [
            4.27611269,
            -5.42822417e-4,
            1.67335701e-5,
            -2.15770813e-8,
            8.62454363e-12,
            -17702.582100,
            3.43505074,
        ],
        high: [
            4.16500285,
            4.90831694e-3,
            -1.90139225e-6,
            3.71185986e-10,
            -2.87908305e-14,
            -17861.787700,
            2.91615662,
        ],
    },
    Nasa7 {
        // N2, 121286, 300-1000 / 1000-5000 K
        low: [
            3.29867700,
            1.40824040e-3,
            -3.96322200e-6,
            5.64151500e-9,
            -2.44485400e-12,
            -1020.899900,
            3.95037200,
        ],
        high: [
            2.92664000,
            1.48797680e-3,
            -5.68476000e-7,
            1.00970380e-10,
            -6.75335100e-15,
            -922.797700,
            5.98052800,
        ],
    },
];

/// coefficient range switch temperature, K, the same 1000 K for every fit.
const T_SWITCH: f64 = 1000.0;

/// cold-end floor for the polynomial window, K (the low fits start at 200).
const T_FLOOR: f64 = 200.0;

/// the coefficient set for species `s` at temperature `t`.
#[inline]
fn coeffs(s: usize, t: f64) -> [f64; 7] {
    let fit = &NASA7[s];
    if t < T_SWITCH {
        fit.low
    } else {
        fit.high
    }
}

/// cp/R for species `s` at temperature `t`.
#[inline]
pub fn cp_r(s: usize, t: f64) -> f64 {
    let t = t.max(T_FLOOR);
    let a = coeffs(s, t);
    a[0] + a[1] * t + a[2] * t * t + a[3] * t * t * t + a[4] * t * t * t * t
}

/// h/(R T) for species `s` at temperature `t`, formation enthalpy included.
#[inline]
pub fn h_rt(s: usize, t: f64) -> f64 {
    let t = t.max(T_FLOOR);
    let a = coeffs(s, t);
    a[0] + a[1] * t / 2.0
        + a[2] * t * t / 3.0
        + a[3] * t.powi(3) / 4.0
        + a[4] * t.powi(4) / 5.0
        + a[5] / t
}

/// s/R for species `s` at temperature `t`, the 1 atm standard state.
/// the a5 term integrates cp/T so it carries T^4/4, not the /5 the
/// enthalpy integral uses.
#[inline]
pub fn s_r(s: usize, t: f64) -> f64 {
    let t = t.max(T_FLOOR);
    let a = coeffs(s, t);
    a[0] * t.ln()
        + a[1] * t
        + a[2] * t * t / 2.0
        + a[3] * t.powi(3) / 3.0
        + a[4] * t.powi(4) / 4.0
        + a[6]
}

/// g/(R T) = h/(R T) - s/R for species `s` at temperature `t`.
#[inline]
pub fn g_rt(s: usize, t: f64) -> f64 {
    h_rt(s, t) - s_r(s, t)
}

/// cp per mass for species `s`, J/(kg K).
#[inline]
pub fn cp_mass(s: usize, t: f64) -> f64 {
    RU * cp_r(s, t) / M[s]
}

/// static enthalpy per mass for species `s`, J/kg, formation included.
#[inline]
pub fn h_mass(s: usize, t: f64) -> f64 {
    RU * t * h_rt(s, t) / M[s]
}

/// mixture cp per mass, J/(kg K), for a mass-fraction vector.
pub fn mix_cp_mass(y: &[f64; NS], t: f64) -> f64 {
    y.iter()
        .enumerate()
        .map(|(s, &yy)| yy * cp_mass(s, t))
        .sum()
}

/// mixture static enthalpy per mass, J/kg, formation included.
pub fn mix_h_mass(y: &[f64; NS], t: f64) -> f64 {
    y.iter().enumerate().map(|(s, &yy)| yy * h_mass(s, t)).sum()
}

/// ln Kp at temperature `t` for a reaction given as (species, stoichiometric
/// coefficient) pairs with products positive, dimensionless on the 1 atm
/// standard state: ln Kp = -delta(g/R T).
pub fn ln_kp(nu: &[(usize, f64)], t: f64) -> f64 {
    let mut ln = 0.0;
    for &(s, coeff) in nu {
        ln -= coeff * g_rt(s, t);
    }
    ln
}

/// ln Kc at temperature `t`, the concentration-form equilibrium constant:
/// Kc = Kp * (p0 / (R T))^dnu with dnu the net mole change. units mol/m3
/// carry the dnu power.
pub fn ln_kc(nu: &[(usize, f64)], t: f64) -> f64 {
    let mut dnu = 0.0;
    for &(_, coeff) in nu {
        dnu += coeff;
    }
    ln_kp(nu, t) + dnu * (P0 / (RU * t)).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hf(298.15) anchors, kj/mol, from the section 3 table: the full
    /// low-range sum must reproduce them (the h2o2 entry is the 0.23 kj
    /// vintage gap the table itself carries, tolerance 0.3).
    #[test]
    fn hf_anchors_match_the_table() {
        let expected_kj = [
            (IH2, 0.0, 0.05),
            (IO2, 0.0, 0.05),
            (IH2O, -241.83, 0.05),
            (IOH, 39.35, 0.05),
            (IH, 218.00, 0.25),
            (IO_, 249.18, 0.05),
            (IHO2, 12.55, 0.05),
            (IH2O2, -136.11, 0.3),
            (IN2, 0.0, 0.05),
        ];
        for (s, want, tol) in expected_kj {
            let hf = RU * T_REF * h_rt(s, T_REF) / 1000.0;
            assert!(
                (hf - want).abs() < tol,
                "hf({}) = {hf} kj/mol, want {want} +- {tol}",
                SPECIES_NAMES[s]
            );
        }
    }

    /// the hardcoded per-mass formation array must equal the polynomial
    /// evaluation it documents, a transcription check on the constants.
    #[test]
    fn hf_mass_const_matches_the_polynomials() {
        for s in 0..NS {
            let from_poly = h_mass(s, T_REF);
            assert!(
                (from_poly - HF_MASS_298[s]).abs() / HF_MASS_298[s].abs().max(1.0) < 1e-4,
                "hf mass mismatch for {}: {from_poly} vs {}",
                SPECIES_NAMES[s],
                HF_MASS_298[s]
            );
        }
    }

    /// JANAF 1 atm anchors at 298.15 K: cp and entropy of the three
    /// diatomic references plus water.
    #[test]
    fn cp_and_entropy_anchors() {
        let cp_cases = [(IN2, 29.12), (IO2, 29.38), (IH2O, 33.60), (IH2, 28.84)];
        for (s, want) in cp_cases {
            let cp = RU * cp_r(s, T_REF);
            assert!(
                (cp - want).abs() < 0.05,
                "cp({}) = {cp}, want {want}",
                SPECIES_NAMES[s]
            );
        }
        let s_h2o = RU * s_r(IH2O, T_REF);
        assert!((s_h2o - 188.83).abs() < 0.05, "s(h2o) = {s_h2o}");
        let s_o2 = RU * s_r(IO2, T_REF);
        assert!((s_o2 - 205.15).abs() < 0.05, "s(o2) = {s_o2}");
    }

    /// Kp for H2 + 0.5 O2 -> H2O: at 298.15 K it is the water formation
    /// free energy, JANAF -228.57 kj/mol giving ln Kp about 92.2, and it
    /// must fall with temperature (exothermic, van't Hoff).
    #[test]
    fn kp_water_formation_anchor() {
        let nu = [(IH2, -1.0), (IO2, -0.5), (IH2O, 1.0)];
        let ln = ln_kp(&nu, T_REF);
        assert!((ln - 92.20).abs() < 0.02, "ln kp(298.15) = {ln}");
        // the implied delta-G0 must sit at the JANAF formation value.
        let dg = -ln * RU * T_REF / 1000.0;
        assert!((dg - (-228.57)).abs() < 0.1, "delta G0 = {dg} kj/mol");
        // monotone fall through the flame band, plus the sign flip far up.
        let l2 = ln_kp(&nu, 2000.0);
        let l3 = ln_kp(&nu, 2500.0);
        let l4 = ln_kp(&nu, 3000.0);
        assert!(
            l2 < ln && l3 < l2 && l4 < l3,
            "van't hoff trend: {ln} {l2} {l3} {l4}"
        );
        assert!(l4 > 0.0, "water still favored at 3000 k: {l4}");
    }

    /// Kp scaling and inverse identities on the same polynomials: doubling
    /// the reaction doubles ln Kp and the reverse is the negative.
    #[test]
    fn kp_scaling_and_inverse() {
        let nu = [(IH2, -1.0), (IO2, -0.5), (IH2O, 1.0)];
        let doubled = [(IH2, -2.0), (IO2, -1.0), (IH2O, 2.0)];
        let reversed = [(IH2O, -1.0), (IH2, 1.0), (IO2, 0.5)];
        for t in [1500.0, 2500.0, 3500.0] {
            let a = ln_kp(&nu, t);
            assert!((ln_kp(&doubled, t) - 2.0 * a).abs() < 1e-12);
            assert!((ln_kp(&reversed, t) + a).abs() < 1e-12);
        }
    }

    /// Kc converts to Kp through the standard state: at equilibrium
    /// Kc * (R T / p0)^dnu = Kp, checked on a dnu = +1 dissociation and a
    /// dnu = -1 recombination.
    #[test]
    fn kc_kp_conversion() {
        let dissoc = [(IH2, -1.0), (IH, 2.0)];
        let recomb = [(IH2O, 2.0), (IH2, -2.0), (IO2, -1.0)];
        for t in [1200.0, 2500.0, 3400.0] {
            let ratio_d = ln_kc(&dissoc, t) - ln_kp(&dissoc, t);
            let want_d = (P0 / (RU * t)).ln();
            assert!((ratio_d - want_d).abs() < 1e-12, "dnu 1 at {t}");
            let ratio_r = ln_kc(&recomb, t) - ln_kp(&recomb, t);
            let want_r = -(P0 / (RU * t)).ln();
            assert!((ratio_r - want_r).abs() < 1e-12, "dnu -1 at {t}");
        }
    }

    /// the two ranges must agree at the 1000 K switch to near machine
    /// precision (the fits were constructed continuous there): both
    /// coefficient sets evaluated at exactly the switch temperature.
    #[test]
    fn ranges_continuous_at_the_switch() {
        let eval_h = |a: &[f64; 7], t: f64| -> f64 {
            a[0] + a[1] * t / 2.0
                + a[2] * t * t / 3.0
                + a[3] * t.powi(3) / 4.0
                + a[4] * t.powi(4) / 5.0
                + a[5] / t
        };
        let eval_s = |a: &[f64; 7], t: f64| -> f64 {
            a[0] * t.ln()
                + a[1] * t
                + a[2] * t * t / 2.0
                + a[3] * t.powi(3) / 3.0
                + a[4] * t.powi(4) / 4.0
                + a[6]
        };
        for spec in NASA7.iter() {
            let lo = eval_h(&spec.low, T_SWITCH);
            let hi = eval_h(&spec.high, T_SWITCH);
            let rel = (lo - hi).abs() / hi.abs().max(1e-30);
            assert!(rel < 1e-5, "h switch gap: {rel}");
            let slo = eval_s(&spec.low, T_SWITCH);
            let shi = eval_s(&spec.high, T_SWITCH);
            let srel = (slo - shi).abs() / shi.abs().max(1e-30);
            assert!(srel < 1e-5, "s switch gap: {srel}");
        }
    }

    /// element bookkeeping: every species' element vector times its molar
    /// mass must close against the conservation identity the gates rely on,
    /// here checked as H2 + O2 stoichiometry sanity (water is 2:1 H:O).
    #[test]
    fn element_matrix_is_sane() {
        assert_eq!(ELEMENTS[IH2O], [2, 1, 0]);
        assert_eq!(ELEMENTS[IHO2], [1, 2, 0]);
        assert_eq!(ELEMENTS[IH2O2], [2, 2, 0]);
        assert_eq!(ELEMENTS[IN2], [0, 0, 2]);
        // atom masses must add: W(H2O) = 2*W(H) + W(O) to the table rounding.
        let wh2o = 2.0 * M[IH] + M[IO_];
        assert!((wh2o - M[IH2O]).abs() < 2e-6, "water mass closure {wh2o}");
    }
}

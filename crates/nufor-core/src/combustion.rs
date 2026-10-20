//! li 2004 h2/o2 mechanism, 19 reactions over 9 species: forward
//! rates, third-body concentrations, fc-only troe falloff, reverse
//! rates through the nasa-7 equilibrium constants, and the species
//! source term omega.
//!
//! rate convention, chemkin: k = A T^n exp(-E / (R T)). li's table i
//! prints cm3, mol, s, kcal, K. the constants below are converted to
//! si on the MOL basis (concentrations mol/m3, masses kg/mol, the
//! same basis `combustion_thermo::M` fixes):
//!
//! - bimolecular rate constants (two concentration factors in the
//!   law): A * 1e-6, m3/mol/s.
//! - termolecular rate constants (three factors): A * 1e-12,
//!   m6/mol2/s. note the factor is 1e-12, not the 1e-6 a kmol-based
//!   table would carry: cm6 = (1e-2 m)^6.
//! - activation energies: E * 4184, J/mol.
//!
//! reaction 5 (H2 + M = H + H + M) is bimolecular-A (law k [H2] [M]);
//! reactions 6, 7, 8 are termolecular-A (laws k [O][O][M], k [O][H][M],
//! k [H][OH][M]).
//!
//! falloff reactions 9 (H + O2 (+M) = HO2 (+M), footnote-c N2 branch)
//! and 15 (H2O2 (+M) = OH + OH (+M)) use the reduced 1-parameter troe
//! form li's table supports: F = Fc constant, so
//!
//! `k(T, [M]) = kinf * (pr / (1 + pr)) * Fc,  pr = k0 [M] / kinf`
//!
//! with the rate law carrying no extra [M]: the third body is inside
//! k through pr. for reaction 9 k0 is termolecular (m6/mol2/s) and
//! kinf bimolecular (m3/mol/s); for reaction 15 k0 is bimolecular
//! (m3/mol/s) and kinf unimolecular (1/s), the standard (+M)
//! chemkin semantics.
//!
//! reverse rates: k_b = k_f / K_c with K_c = K_p (p0 / (R T))^dnu
//! from the delta-G0 of the same nasa-7 polynomials (the standard
//! chemkin procedure; the mechanism's published validation runs were
//! built on exactly this pairing, so the fits are used unpatched).
//!
//! the oh thermo carries the documented vintage bias (Hf = 39.35
//! kj/mol, about 0.5 kj/mol above the modern 8.91 kcal/mol in li's
//! table ii), so equilibrium X_OH reads low against a modern cea run;
//! that bias is inside the gate tolerances, never patched here.

use crate::combustion_thermo::{
    h_mass, ln_kc, IH, IH2, IH2O, IH2O2, IHO2, IO2, IOH, IO_, M, NS, RU,
};

/// number of reactions in the li 2004 set.
pub const NR: usize = 19;

/// forward arrhenius triples (A, n, E) in si units, one per reaction;
/// rows 8 and 14 are falloff (see `FALLOFF`) and carry zeros here.
/// reactions 14 and 19 are duplicate pairs whose two expressions sum;
/// the second expression sits in `R14_SECOND` / `R19_SECOND`.
pub const RATE_TRIPLES: [(f64, f64, f64); NR] = [
    (3.55e15 * 1e-6, -0.41, 16.6 * 4184.0),   // 1: H + O2 = O + OH
    (5.08e4 * 1e-6, 2.67, 6.29 * 4184.0),     // 2: O + H2 = H + OH
    (2.16e8 * 1e-6, 1.51, 3.43 * 4184.0),     // 3: H2 + OH = H2O + H
    (2.97e6 * 1e-6, 2.02, 13.4 * 4184.0),     // 4: O + H2O = OH + OH
    (4.58e19 * 1e-6, -1.40, 104.38 * 4184.0), // 5: H2 + M = H + H + M, bimolecular-A
    (6.16e15 * 1e-12, -0.50, 0.0),            // 6: O + O + M = O2 + M, termolecular-A
    (4.71e18 * 1e-12, -1.00, 0.0),            // 7: O + H + M = OH + M, termolecular-A
    (3.8e22 * 1e-12, -2.00, 0.0),             // 8: H + OH + M = H2O + M, termolecular-A
    (0.0, 0.0, 0.0),                          // 9: H + O2 (+M) = HO2 (+M), falloff
    (1.66e13 * 1e-6, 0.00, 0.82 * 4184.0),    // 10: HO2 + H = H2 + O2
    (7.08e13 * 1e-6, 0.00, 0.30 * 4184.0),    // 11: HO2 + H = OH + OH
    (3.25e13 * 1e-6, 0.00, 0.0),              // 12: HO2 + O = OH + O2
    (2.89e13 * 1e-6, 0.00, -0.50 * 4184.0),   // 13: HO2 + OH = H2O + O2
    (4.20e14 * 1e-6, 0.00, 11.98 * 4184.0),   // 14a: HO2 + HO2 = H2O2 + O2
    (0.0, 0.0, 0.0),                          // 15: H2O2 (+M) = OH + OH (+M), falloff
    (2.41e13 * 1e-6, 0.00, 3.97 * 4184.0),    // 16: H2O2 + H = H2O + OH
    (4.82e13 * 1e-6, 0.00, 7.95 * 4184.0),    // 17: H2O2 + H = H2 + HO2
    (9.55e6 * 1e-6, 2.00, 3.97 * 4184.0),     // 18: H2O2 + O = OH + HO2
    (1.00e12 * 1e-6, 0.00, 0.0),              // 19a: H2O2 + OH = H2O + HO2
];

/// arrhenius triple, si units: (pre-exponential, temperature exponent, activation j/mol).
pub type ArrheniusTriple = (f64, f64, f64);

/// reaction 14's second (duplicate) expression, si units.
pub const R14_SECOND: ArrheniusTriple = (1.30e11 * 1e-6, 0.00, -1.63 * 4184.0);

/// reaction 19's second (duplicate) expression, si units.
pub const R19_SECOND: ArrheniusTriple = (5.8e14 * 1e-6, 0.00, 9.56 * 4184.0);

/// troe falloff data, si: (reaction index, k0 triple, kinf triple, Fc).
/// reaction 9: k0 termolecular (footnote-c N2 branch), kinf bimolecular.
/// reaction 15: k0 bimolecular, kinf unimolecular (1/s).
pub const FALLOFF: [(usize, ArrheniusTriple, ArrheniusTriple, f64); 2] = [
    (
        8,
        (6.37e20 * 1e-12, -1.72, 0.52 * 4184.0),
        (1.48e12 * 1e-6, 0.60, 0.0),
        0.8,
    ),
    (
        14,
        (1.20e17 * 1e-6, 0.00, 45.5 * 4184.0),
        (2.95e14, 0.00, 48.4 * 4184.0),
        0.5,
    ),
];

/// dense stoichiometry, products positive; the third body is never a
/// row entry (it cancels in every reaction of this mechanism).
pub const NU_TABLE: [[f64; NS]; NR] = build_nu_table();

const fn build_nu_table() -> [[f64; NS]; NR] {
    let mut table = [[0.0f64; NS]; NR];
    // (reaction index, species index, coefficient)
    let rows: [(usize, usize, f64); 64] = [
        (0, IH, -1.0),
        (0, IO2, -1.0),
        (0, IO_, 1.0),
        (0, IOH, 1.0),
        (1, IO_, -1.0),
        (1, IH2, -1.0),
        (1, IH, 1.0),
        (1, IOH, 1.0),
        (2, IH2, -1.0),
        (2, IOH, -1.0),
        (2, IH2O, 1.0),
        (2, IH, 1.0),
        (3, IO_, -1.0),
        (3, IH2O, -1.0),
        (3, IOH, 2.0),
        (4, IH2, -1.0),
        (4, IH, 2.0),
        (5, IO_, -2.0),
        (5, IO2, 1.0),
        (6, IO_, -1.0),
        (6, IH, -1.0),
        (6, IOH, 1.0),
        (7, IH, -1.0),
        (7, IOH, -1.0),
        (7, IH2O, 1.0),
        (8, IH, -1.0),
        (8, IO2, -1.0),
        (8, IHO2, 1.0),
        (9, IHO2, -1.0),
        (9, IH, -1.0),
        (9, IH2, 1.0),
        (9, IO2, 1.0),
        (10, IHO2, -1.0),
        (10, IH, -1.0),
        (10, IOH, 2.0),
        (11, IHO2, -1.0),
        (11, IO_, -1.0),
        (11, IOH, 1.0),
        (11, IO2, 1.0),
        (12, IHO2, -1.0),
        (12, IOH, -1.0),
        (12, IH2O, 1.0),
        (12, IO2, 1.0),
        (13, IHO2, -2.0),
        (13, IH2O2, 1.0),
        (13, IO2, 1.0),
        (14, IH2O2, -1.0),
        (14, IOH, 2.0),
        (15, IH2O2, -1.0),
        (15, IH, -1.0),
        (15, IH2O, 1.0),
        (15, IOH, 1.0),
        (16, IH2O2, -1.0),
        (16, IH, -1.0),
        (16, IH2, 1.0),
        (16, IHO2, 1.0),
        (17, IH2O2, -1.0),
        (17, IO_, -1.0),
        (17, IOH, 1.0),
        (17, IHO2, 1.0),
        (18, IH2O2, -1.0),
        (18, IOH, -1.0),
        (18, IH2O, 1.0),
        (18, IHO2, 1.0),
    ];
    let mut i = 0;
    while i < rows.len() {
        table[rows[i].0][rows[i].1] = rows[i].2;
        i += 1;
    }
    table
}

/// sparse stoichiometry for the equilibrium-constant calls, at most 4
/// nonzero entries per reaction, padded with (0, 0.0) which contributes
/// nothing to ln Kp or dnu.
pub const NU_PAIRS: [[(usize, f64); 4]; NR] = build_nu_pairs();

const fn build_nu_pairs() -> [[(usize, f64); 4]; NR] {
    let mut out = [[(0usize, 0.0f64); 4]; NR];
    let mut r = 0;
    while r < NR {
        let mut s = 0;
        let mut j = 0;
        while s < NS {
            let nu = NU_TABLE[r][s];
            if nu != 0.0 {
                out[r][j] = (s, nu);
                j += 1;
            }
            s += 1;
        }
        r += 1;
    }
    out
}

/// whether reaction `r` is a troe falloff reaction (9 or 15).
#[inline]
pub fn is_falloff(r: usize) -> bool {
    matches!(r, 8 | 14)
}

/// whether reaction `r` takes a third body in its rate law (the plain
/// third-body reactions; the falloff pair takes [M] through pr instead).
#[inline]
pub fn needs_third_body(r: usize) -> bool {
    matches!(r, 4..=7)
}

/// third-body efficiencies from li's footnotes a, b, c, f; default 1.0
/// (N2 carries the default). Ar and He efficiencies exist in the source
/// but those species are absent from the 9-species set.
#[inline]
pub fn third_body_eff(r: usize, s: usize) -> f64 {
    match (r, s) {
        (4, IH2O) | (5, IH2O) | (6, IH2O) => 12.0, // footnote a
        (4, IH2) | (5, IH2) | (6, IH2) => 2.5,
        (7, IH2O) => 12.0, // footnote b
        (7, IH2) => 2.5,
        (8, IH2O) => 11.0, // footnote c, N2 branch
        (8, IH2) => 2.0,
        (8, IO2) => 0.78,
        (14, IH2O) => 12.0, // footnote f
        (14, IH2) => 2.5,
        _ => 1.0,
    }
}

/// arrhenius evaluation k = A T^n exp(-E / (R T)); a negative E is legal
/// (chemkin semantics, reaction 13 and 14b carry one).
#[inline]
pub fn arrhenius(a: f64, n: f64, e: f64, t: f64) -> f64 {
    a * t.powf(n) * (-e / (RU * t)).exp()
}

/// forward rate constant for a non-falloff reaction `r` at temperature
/// `t`, si units. falloff rows return 0.0: their constant is undefined
/// without the third-body concentration, see `k_falloff`.
pub fn k_f(r: usize, t: f64) -> f64 {
    if is_falloff(r) {
        return 0.0;
    }
    let (a, n, e) = RATE_TRIPLES[r];
    let mut k = arrhenius(a, n, e, t);
    if r == 13 {
        let (a2, n2, e2) = R14_SECOND;
        k += arrhenius(a2, n2, e2, t);
    }
    if r == 18 {
        let (a2, n2, e2) = R19_SECOND;
        k += arrhenius(a2, n2, e2, t);
    }
    k
}

/// concentration-form equilibrium constant K_c(r, t) from the same
/// nasa-7 fits the energy closure uses; units carry the dnu power
/// (mol/m3 per unit dnu).
pub fn k_c(r: usize, t: f64) -> f64 {
    ln_kc(&NU_PAIRS[r], t).exp()
}

/// backward rate constant for a non-falloff reaction, k_b = k_f / K_c.
pub fn k_b(r: usize, t: f64) -> f64 {
    k_f(r, t) / k_c(r, t)
}

/// the fc-only troe falloff constant at temperature `t` with third-body
/// concentration `cm` (mol/m3, efficiency-weighted).
pub fn k_falloff(r: usize, t: f64, cm: f64) -> f64 {
    let (_, (a0, n0, e0), (ainf, ninf, einf), fc) = FALLOFF[if r == 8 { 0 } else { 1 }];
    let k0 = arrhenius(a0, n0, e0, t);
    let kinf = arrhenius(ainf, ninf, einf, t);
    if kinf == 0.0 {
        return 0.0;
    }
    let pr = k0 * cm / kinf;
    kinf * (pr / (1.0 + pr)) * fc
}

/// number densities, mol/m3, from (rho, Y).
pub fn conc_species(rho: f64, y: &[f64; NS]) -> [f64; NS] {
    let n: [f64; NS] = core::array::from_fn(|s| rho * y[s] / M[s]);
    n
}

/// efficiency-weighted third-body concentration for reaction `r`,
/// mol/m3. plain third-body reactions and both falloff reactions use it
/// (the falloff pair through the reduced pressure pr).
pub fn conc_m(r: usize, n: &[f64; NS]) -> f64 {
    let m: f64 = (0..NS).map(|s| third_body_eff(r, s) * n[s]).sum();
    m
}

/// net rate of reaction `r`, mol/m3/s, forward minus reverse. the
/// reverse constant is k_f / K_c everywhere, the falloff pair included
/// (the pressure blend commutes with reversibility, so detailed balance
/// holds at the blended constant).
fn reaction_rate(r: usize, n: &[f64; NS], t: f64) -> f64 {
    let kc = k_c(r, t);
    let kf = if is_falloff(r) {
        k_falloff(r, t, conc_m(r, n))
    } else {
        k_f(r, t)
    };
    if kf == 0.0 {
        return 0.0;
    }
    let mut fwd = kf;
    let mut rev = kf / kc;
    for (s, &nu) in NU_TABLE[r].iter().enumerate() {
        if nu < 0.0 {
            fwd *= n[s].powf(-nu);
        } else if nu > 0.0 {
            rev *= n[s].powf(nu);
        }
    }
    if needs_third_body(r) {
        let cm = conc_m(r, n);
        fwd *= cm;
        rev *= cm;
    }
    fwd - rev
}

/// one cell's reaction state: (rho, Y, T) gives the per-species mass
/// production d(Y_i)/dt, 1/s, and the chemical energy deposition
/// d(e_chem)/dt, J/(kg s) — the mixture static-enthalpy rate the
/// composition change carries at fixed T, so a positive value means the
/// composition is climbing the enthalpy ladder and draining the thermal
/// pool (the march keeps total energy in e and recovers T through the
/// closure, the same contract `species_chem::omega` serves).
pub fn omega(rho: f64, y: &[f64; NS], t: f64) -> ([f64; NS], f64) {
    let n = conc_species(rho, y);
    let mut dydt = [0.0f64; NS];
    for (r, row) in NU_TABLE.iter().enumerate() {
        let rate = reaction_rate(r, &n, t);
        if rate == 0.0 {
            continue;
        }
        for s in 0..NS {
            let nu = row[s];
            if nu != 0.0 {
                dydt[s] += nu * rate * M[s] / rho;
            }
        }
    }
    let de: f64 = dydt.iter().enumerate().map(|(s, d)| d * h_mass(s, t)).sum();
    (dydt, de)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::combustion_thermo::SPECIES_NAMES;

    /// forward anchors at two temperatures, hand-computed from li's
    /// table i values (A, n, E printed there) with the si conversions:
    /// k = A T^n exp(-E/(R T)), R = 8.314462. 1500 K and 2500 K, both
    /// inside the flame band. hand arithmetic carries about 4 digits,
    /// so the tolerance is 1 percent.
    #[test]
    fn forward_rate_anchors() {
        let cases: [(usize, f64, f64); 8] = [
            // r1 H + O2 = O + OH: 3.55e9 * 1500^-0.41 * exp(-69454.4/RT)
            (0, 1500.0, 6.7524e5),
            (0, 2500.0, 5.0808e6),
            // r2 O + H2 = H + OH: 5.08e-2 * T^2.67 * exp(-26317.4/RT)
            (1, 1500.0, 1.8603e6),
            (1, 2500.0, 1.6924e7),
            // r3 H2 + OH = H2O + H: 2.16e2 * T^1.51 * exp(-14351.1/RT)
            (2, 1500.0, 4.2718e6),
            (2, 2500.0, 1.4639e7),
            // r5 H2 + M = 2H + M, bimolecular-A: 4.58e13 * T^-1.4 * exp(-436726/RT)
            (4, 1500.0, 1.0150e-6),
            (4, 2500.0, 6.0123e-1),
        ];
        for (r, t, want) in cases {
            let k = k_f(r, t);
            let rel = (k - want).abs() / want;
            assert!(rel < 0.01, "k_f({r}) at {t} = {k}, want {want}, rel {rel}");
        }
    }

    /// termolecular-unit anchors: reactions 6 and 8 carry three
    /// concentration factors, so A converts with 1e-12 to m6/mol2/s.
    /// r6 at 1500 K: 6.16e3 * 1500^-0.5 = 159.0; at 2500 K: 123.2.
    /// r8 at 1500 K: 3.8e10 * 1500^-2 = 1.689e4; at 2500 K: 6.08e3.
    #[test]
    fn termolecular_anchors() {
        let cases: [(usize, f64, f64); 4] = [
            (5, 1500.0, 159.0),
            (5, 2500.0, 123.2),
            (7, 1500.0, 1.689e4),
            (7, 2500.0, 6.08e3),
        ];
        for (r, t, want) in cases {
            let k = k_f(r, t);
            let rel = (k - want).abs() / want;
            assert!(rel < 0.01, "k_f({r}) at {t} = {k}, want {want}, rel {rel}");
        }
    }

    /// duplicate-pair anchors, the two expressions summing:
    /// r14 at 1500 K: 4.2e8 exp(-50124.3/RT) + 1.3e5 exp(6819.9/RT) =
    /// 7.78e6; r19 at 2500 K: 1e6 + 5.8e8 exp(-39999/RT) = 8.57e7.
    #[test]
    fn duplicate_rates_sum() {
        let cases: [(usize, f64, f64); 2] = [(13, 1500.0, 7.7720e6), (18, 2500.0, 8.5666e7)];
        for (r, t, want) in cases {
            let k = k_f(r, t);
            let rel = (k - want).abs() / want;
            assert!(rel < 0.01, "k_f({r}) at {t} = {k}, want {want}, rel {rel}");
        }
    }

    /// falloff anchors at 1500 K with [M] = 1 mol/m3, hand-computed:
    /// r9 k0 = 6.37e8 * 1500^-1.72 * exp(-2175.7/RT) = 1848 m6/mol2/s,
    /// kinf = 1.48e6 * 1500^0.6 = 1.191e8 m3/mol/s, pr = 1.552e-5,
    /// k = kinf (pr/(1+pr)) 0.8 = 1479 m3/mol/s.
    /// r15 k0 = 1.2e11 exp(-190372/RT) = 2.804e4 m3/mol/s, kinf =
    /// 2.95e14 exp(-202505.6/RT) = 2.639e7 1/s, pr = 1.066e-3,
    /// k = kinf (pr/(1+pr)) 0.5 = 1.405e4 1/s.
    #[test]
    fn falloff_anchors() {
        let k9 = k_falloff(8, 1500.0, 1.0);
        assert!((k9 - 1479.0).abs() / 1479.0 < 0.01, "r9 blended k = {k9}");
        let k15 = k_falloff(14, 1500.0, 1.0);
        assert!(
            (k15 - 1.405e4).abs() / 1.405e4 < 0.01,
            "r15 blended k = {k15}"
        );
        // the blend saturates: at [M] = 1e6 mol/m3 pr ~ 1e6*... and k -> Fc*kinf.
        let kinf9 = arrhenius(1.48e12 * 1e-6, 0.60, 0.0, 1500.0);
        let k9_hi = k_falloff(8, 1500.0, 1e10);
        assert!(k9_hi > 0.99 * 0.8 * kinf9, "r9 high-p limit: {k9_hi}");
        // and vanishes with the third body.
        assert_eq!(k_falloff(8, 1500.0, 0.0), 0.0);
    }

    /// element conservation: every stoichiometry row weighted by molar
    /// masses must close, and the sparse pairs must match the dense table.
    #[test]
    fn stoichiometry_conserves_elements() {
        let el = crate::combustion_thermo::ELEMENTS;
        for (r, row) in NU_TABLE.iter().enumerate() {
            let mut h = 0.0;
            let mut o = 0.0;
            let mut n = 0.0;
            for s in 0..NS {
                h += row[s] * el[s][0] as f64;
                o += row[s] * el[s][1] as f64;
                n += row[s] * el[s][2] as f64;
            }
            assert!(h.abs() < 1e-12, "reaction {r} breaks H balance");
            assert!(o.abs() < 1e-12, "reaction {r} breaks O balance");
            assert!(n.abs() < 1e-12, "reaction {r} breaks N balance");
            // mass closure through molar masses.
            let mut mass = 0.0;
            for s in 0..NS {
                mass += row[s] * M[s];
            }
            assert!(mass.abs() < 1e-15, "reaction {r} breaks mass: {mass}");
        }
        // sparse pairs equal the dense rows.
        for (r, dense) in NU_TABLE.iter().enumerate() {
            let mut dense = *dense;
            for &(s, nu) in NU_PAIRS[r].iter() {
                if s > 0 || nu != 0.0 {
                    assert!((dense[s] - nu).abs() < 1e-15, "pair mismatch r{s}");
                    dense[s] = 0.0;
                }
            }
            assert!(
                dense.iter().all(|&v| v == 0.0),
                "pairs miss entries in r{r}"
            );
        }
    }

    /// omega at a burning state: the mixture moves toward water, mass is
    /// conserved, and the thermal pool heats (negative deposition = the
    /// composition is climbing DOWN the enthalpy ladder = exothermic).
    /// the state is the stoichiometric mix with an oh seed, renormalized
    /// to unit mass fractions (a raw unburned mix carries no radicals and
    /// nearly freezes; a rich mix past radical equilibrium reads
    /// endothermic by design).
    #[test]
    fn omega_moves_toward_water_and_conserves_mass() {
        let t = 2500.0f64;
        let rho = 1.0f64;
        let y_h2 = 2.01588 / (2.01588 + 0.5 * 31.9988);
        let mut y = [0.0f64; NS];
        y[IH2] = y_h2;
        y[IO2] = 1.0 - y_h2;
        y[IOH] = 0.005;
        let sum: f64 = y.iter().sum();
        for v in y.iter_mut() {
            *v /= sum;
        }
        let (dydt, de) = omega(rho, &y, t);
        let largest = dydt.iter().fold(0.0f64, |a, d| a.max(d.abs()));
        let sum: f64 = dydt.iter().sum();
        assert!(sum.abs() < 1e-10 * largest, "mass not conserved: {sum}");
        assert!(dydt[IH2] < 0.0, "H2 must burn");
        assert!(dydt[IO2] < 0.0, "O2 must burn");
        assert!(dydt[IH2O] > 0.0, "H2O must grow");
        assert!(de < 0.0, "exothermic deposition: {de}");
        for (s, d) in dydt.iter().enumerate() {
            assert!(d.is_finite(), "{} diverged: {}", SPECIES_NAMES[s], d);
        }
    }

    /// detailed balance: a state holding only H2 and H leaves exactly one
    /// reaction live (r5, H2 + M = H + H + M; every other reaction needs
    /// an O-bearing or HO2/H2O2 species). at its Kc-derived equilibrium
    /// composition n(H2) = n(H)^2 / Kc the net rate must vanish to
    /// machine precision, which pins the k_b = k_f / Kc wiring.
    #[test]
    fn equilibrium_nulls_the_isolated_h2_h_pair() {
        let t = 2500.0f64;
        let rho = 1.0f64;
        let kc5 = k_c(4, t);
        let n_h = 0.02f64;
        let n_h2 = n_h * n_h / kc5;
        let mut y = [0.0f64; NS];
        y[IH2] = n_h2 * M[IH2] / rho;
        y[IH] = n_h * M[IH] / rho;
        let (dydt, _) = omega(rho, &y, t);
        // the one-way flux scales the residual: near equilibrium the
        // net is a difference of two nearly equal large fluxes.
        let n = conc_species(rho, &y);
        let flux = k_f(4, t) * n[IH2] * conc_m(4, &n);
        let net = dydt[IH2].abs() * rho / M[IH2];
        assert!(
            (net / flux).abs() < 1e-12,
            "r5 not nulled at equilibrium: net/flux = {}",
            net / flux
        );
        assert!(dydt[IHO2].abs() < 1e-30, "ho2 must stay frozen");
        assert!(dydt[IO2].abs() < 1e-30, "o2 must stay frozen");
    }

    // ---- cea gates -------------------------------------------------------
    //
    // gate 2 (cea_out/h2o2_tp_*.out): equilibrium composition at assigned
    // (t, p). omega's fixed point is the detailed-balance state: every
    // reaction's forward and reverse fluxes equalize there, and through
    // k_b = k_f / K_c that is prod(n^nu) = Kc with Kc from the same
    // nasa-7 polynomials the reverse rates use. the gate solves those
    // algebraic conditions by newton on the element potentials (the
    // standard equilibrium-constant method) and then verifies omega
    // really is nulled at the solution, which ties the composition to
    // the rate machinery this module actually ships.
    //
    // integrating omega in time to the same state stalls: the peroxide
    // pool is quasi-steady at 2500 k, 70 bar (the r14-r19 balance sits
    // at 1.6x for hundreds of thousands of explicit steps while the
    // mixture composition has already stopped moving at the 1e-4 level;
    // measured while building this gate). the algebraic solve is the
    // convergent route to the same fixed point.
    //
    // the oh column reads ~7-9% low against cea 2024 across every state:
    // the documented thermo-vintage bias (h2o2.md section 4, the rus78
    // fit carries hf(oh) = 39.35 kj/mol vs the modern 8.91 kcal/mol).
    // that bias is IN tolerance here, never patched.

    /// element-potential newton for the detailed-balance composition at
    /// (t, p): unknowns lam_h, lam_o, ln_ntot (nitrogen absent from the
    /// gate states), equations: o/h element ratio, ideal-gas total moles,
    /// and the 6 stoichiometrically independent reaction balances
    /// (reactions 1, 2, 3, 5, 9, 15, one-indexed). returns mole
    /// fractions and the max residual.
    fn equilibrium_tp(t: f64, p: f64, o_per_h: f64) -> ([f64; NS], f64) {
        use crate::combustion_thermo::{g_rt, ln_kc, ELEMENTS};
        // the 6 independent rows as (species, nu) pairs.
        let nu_rows: [[(usize, f64); 4]; 6] = [
            NU_PAIRS[0],
            NU_PAIRS[1],
            NU_PAIRS[2],
            NU_PAIRS[4],
            NU_PAIRS[8],
            NU_PAIRS[14],
        ];
        // active species: all but N2 (the gate mixtures carry no nitrogen).
        let act: [usize; 8] = [IH2, IO2, IH2O, IOH, IH, IO_, IHO2, IH2O2];
        let g: [f64; NS] = core::array::from_fn(|s| g_rt(s, t));
        // v = (lam_h, lam_o, ln ntot); n_s = ntot exp(x_s)/Z with
        // x_s = -g_s - lam_h a_hs - lam_o a_os.
        let pot = |v: [f64; 3], s: usize| -> f64 {
            -g[s] - v[0] * ELEMENTS[s][0] as f64 - v[1] * ELEMENTS[s][1] as f64
        };
        let comp = |v: [f64; 3]| -> ([f64; NS], f64) {
            // returns (mole fractions, z)
            let mut raw = [0.0f64; NS];
            for &s in act.iter() {
                let e = pot(v, s);
                raw[s] = e.clamp(-690.0, 690.0).exp();
            }
            let z: f64 = act.iter().map(|&s| raw[s]).sum();
            let mut x = [0.0f64; NS];
            for &s in act.iter() {
                x[s] = raw[s] / z;
            }
            (x, z)
        };
        let resid = |v: [f64; 3]| -> [f64; 8] {
            let (x, z) = comp(v);
            let mut out = [0.0f64; 8];
            // o/h ratio (per kg: (ntot/z) cancels in the ratio).
            let mut h_tot = 0.0;
            let mut o_tot = 0.0;
            for &s in act.iter() {
                h_tot += ELEMENTS[s][0] as f64 * x[s];
                o_tot += ELEMENTS[s][1] as f64 * x[s];
            }
            out[0] = (o_tot / h_tot).ln() - o_per_h.ln();
            // ideal-gas mole sum at fixed (t, p): ntot * R T = p.
            out[1] = v[2] + (RU * t).ln() - p.ln();
            // reaction balances: sum nu (ln ntot + pot_s - ln z) = ln Kc.
            for (i, row) in nu_rows.iter().enumerate() {
                let mut e = 0.0f64;
                for &(s, nu) in row.iter() {
                    if nu != 0.0 {
                        e += nu * (v[2] + pot(v, s) - z.ln());
                    }
                }
                out[2 + i] = e - ln_kc(row, t);
            }
            out
        };
        // newton with numeric jacobian and backtracking.
        let mut v = [0.0f64; 3];
        v[2] = (p / (RU * t)).ln().max(1.0);
        let mut f = resid(v);
        for _it in 0..200 {
            let e0 = f.iter().fold(0.0f64, |a, x| a.max(x.abs()));
            if e0 < 1e-13 {
                break;
            }
            let mut jac = [[0.0f64; 3]; 8];
            for j in 0..3 {
                let mut vp = v;
                let h = 1e-7 * (1.0 + v[j].abs());
                vp[j] += h;
                let fp = resid(vp);
                for i in 0..8 {
                    jac[i][j] = (fp[i] - f[i]) / h;
                }
            }
            // normal equations (3x3): solve (j^t j) dv = -j^t f.
            let mut a = [[0.0f64; 4]; 3];
            for i in 0..3 {
                for k in 0..3 {
                    a[i][k] = (0..8).map(|r| jac[r][i] * jac[r][k]).sum();
                }
                a[i][3] = -(0..8).map(|r| jac[r][i] * f[r]).sum::<f64>();
            }
            for c in 0..3 {
                let piv = (c..3).fold(c, |best, r| {
                    if a[r][c].abs() > a[best][c].abs() {
                        r
                    } else {
                        best
                    }
                });
                a.swap(c, piv);
                if a[c][c].abs() < 1e-300 {
                    break;
                }
                for r in 0..3 {
                    if r != c {
                        let fac = a[r][c] / a[c][c];
                        let src = a[c];
                        for (val, s) in a[r].iter_mut().zip(src.iter()).skip(c) {
                            *val -= fac * *s;
                        }
                    }
                }
            }
            let dv = [a[0][3] / a[0][0], a[1][3] / a[1][1], a[2][3] / a[2][2]];
            // backtracking on the max residual.
            let mut dmp = 1.0f64;
            let mut v_new = v;
            while dmp > 1e-10 {
                let vt = [
                    (v[0] + dmp * dv[0]).clamp(-60.0, 60.0),
                    (v[1] + dmp * dv[1]).clamp(-60.0, 60.0),
                    (v[2] + dmp * dv[2]).clamp(-40.0, 40.0),
                ];
                let ft = resid(vt);
                let et = ft.iter().fold(0.0f64, |a, x| a.max(x.abs()));
                if et < e0 {
                    v_new = vt;
                    f = ft;
                    break;
                }
                dmp *= 0.5;
            }
            if v_new == v {
                break;
            }
            v = v_new;
        }
        let (x, _z) = comp(v);
        let e_final = resid(v).iter().fold(0.0f64, |a, x| a.max(x.abs()));
        (x, e_final)
    }

    /// gate 2: the fixed-(t, p) equilibrium compositions against the cea
    /// runs (h2o2_tp_t*.out, stoichiometric h2/o2, no nitrogen). parsing
    /// is python-done once; the values here are transcribed from those
    /// files with the filename in the comment. tolerances: max(0.03
    /// absolute, 15% relative), oh max(0.05 absolute, 25% relative).
    #[test]
    fn cea_tp_equilibrium_composition() {
        // (t, p_bar, cea x per species in index order h2 o2 h2o oh h o)
        // from cea_out/h2o2_tp_*.out (trace=1e-10 runs; ho2/h2o2 excluded:
        // they sit at 1e-7..1e-4, below the gate's reportable floor).
        let states: [(f64, f64, [f64; 6]); 12] = [
            // h2o2_tp_t2500_p1.out
            (
                2500.0,
                1.0,
                [0.043124, 0.015668, 0.908999, 0.025159, 0.00523, 0.001813],
            ),
            // h2o2_tp_t2500_p70.out
            (
                2500.0,
                70.0,
                [0.011052, 0.003946, 0.978174, 0.006392, 0.0003165, 0.0001087],
            ),
            // h2o2_tp_t3000_p1.out
            (
                3000.0,
                1.0,
                [0.134709, 0.045062, 0.639058, 0.099068, 0.058046, 0.02402],
            ),
            // h2o2_tp_t3000_p70.out
            (
                3000.0,
                70.0,
                [0.041838, 0.013501, 0.90894, 0.03022, 0.003866, 0.001571],
            ),
            // h2o2_tp_t3400_p1.out
            (
                3400.0,
                1.0,
                [0.185881, 0.058575, 0.306289, 0.155286, 0.20194, 0.091972],
            ),
            // h2o2_tp_t3400_p70.out
            (
                3400.0,
                70.0,
                [0.086452, 0.025943, 0.79318, 0.070478, 0.01646, 0.007316],
            ),
            // h2o2_tp_t3676.53_p1.out
            (
                3676.53,
                1.0,
                [0.163107, 0.04912, 0.125845, 0.145354, 0.349605, 0.166923],
            ),
            // h2o2_tp_t3676.53_p70.out
            (
                3676.53,
                70.0,
                [0.123872, 0.03546, 0.679398, 0.107625, 0.036415, 0.016951],
            ),
            // h2o2_tp_t4000_p1.out
            (
                4000.0,
                1.0,
                [0.099092, 0.027835, 0.029531, 0.092888, 0.502458, 0.248176],
            ),
            // h2o2_tp_t4000_p70.out
            (
                4000.0,
                70.0,
                [0.165129, 0.044805, 0.522371, 0.152132, 0.077525, 0.037634],
            ),
            // h2o2_tp_t4500_p1.out
            (
                4500.0,
                1.0,
                [0.031824, 0.00776, 0.002154, 0.03088, 0.617679, 0.309698],
            ),
            // h2o2_tp_t4500_p70.out
            (
                4500.0,
                70.0,
                [0.19747, 0.049268, 0.281787, 0.193821, 0.183902, 0.093269],
            ),
        ];
        let o_per_h = 0.5; // stoichiometric: o atoms per h atom = (2*0.5)/2
        let mut checked = 0;
        for (t, p_bar, cea) in states.iter() {
            let (x, res) = equilibrium_tp(*t, p_bar * 1e5, o_per_h);
            assert!(
                res < 1e-10,
                "tp solver did not converge at {t} k: res {res}"
            );
            // omega must be nulled at this composition (ties the gate to
            // the shipped rate machinery): scale by the largest source.
            let ntot = p_bar * 1e5 / (RU * t);
            let mmix: f64 = (0..NS).map(|s| x[s] * M[s]).sum();
            let rho = ntot * mmix;
            let (dydt, _) = omega(rho, &mass_from_x(&x), *t);
            let worst_src = dydt.iter().fold(0.0f64, |a, d| a.max(d.abs()));
            assert!(
                worst_src < 1e-4 * rho,
                "omega not nulled at tp equilibrium {t} k: {worst_src} vs rho {rho}"
            );
            for (s, &sp) in [IH2, IO2, IH2O, IOH, IH, IO_].iter().enumerate() {
                let want = cea[s];
                let tol = if sp == IOH {
                    0.05_f64.max(0.25 * want) // oh vintage bias
                } else {
                    0.03_f64.max(0.15 * want)
                };
                assert!(
                    (x[sp] - want).abs() < tol,
                    "{} at {t} k, {p_bar} bar: nufor {:.5} vs cea {want:.5} (tol {tol:.5})",
                    SPECIES_NAMES[sp],
                    x[sp]
                );
                checked += 1;
            }
        }
        assert_eq!(checked, 72);
    }

    /// gate 1: adiabatic flame temperature, the constant-enthalpy
    /// hp states (cea_out/h2o2_hp_of*_70bar.out). solve t so that
    /// h_mix(x_eq(t), t) = h_reactants (h2 + o2 at 298.15), with x_eq(t)
    /// from the same detailed-balance solver the tp gate uses. cea
    /// tolerance 15 k; the o/f sweep trend must peak near o/f 8.
    #[test]
    fn cea_hp_flame_temperatures() {
        use crate::combustion_thermo::h_mass;
        // (o2 moles per 1 h2, cea t k, cea file)
        let cases: [(f64, f64, &str); 8] = [
            (0.2520, 3137.92, "h2o2_hp_of4_70bar.out"),
            (0.3465, 3521.33, "h2o2_hp_of5.5_70bar.out"),
            (0.3780, 3590.55, "h2o2_hp_of6_70bar.out"),
            (0.4410, 3663.82, "h2o2_hp_of7_70bar.out"),
            (0.5040, 3676.04, "h2o2_hp_of8_70bar.out"),
            (0.6300, 3614.12, "h2o2_hp_of10_70bar.out"),
            (0.7560, 3515.46, "h2o2_hp_of12_70bar.out"),
            (0.5000, 3676.53, "stoich (of2 file)"),
        ];
        let p = 70.0e5;
        let mut temps = Vec::new();
        for (o2_mol, t_cea, file) in cases.iter() {
            let m_tot = M[IH2] + o2_mol * M[IO2];
            let y_h2 = M[IH2] / m_tot;
            let y_o2 = o2_mol * M[IO2] / m_tot;
            let h_react = y_h2 * h_mass(IH2, 298.15) + y_o2 * h_mass(IO2, 298.15);
            let o_per_h = (2.0 * o2_mol) / 2.0;
            // secant on h_eq(t) - h_react.
            let hh = |t: f64| -> f64 {
                let (x, res) = equilibrium_tp(t, p, o_per_h);
                assert!(res < 1e-10, "tp solver failed in flame solve at {t}: {res}");
                let mmix: f64 = (0..NS).map(|s| x[s] * M[s]).sum();
                let h: f64 = (0..NS).map(|s| x[s] * M[s] / mmix * h_mass(s, t)).sum();
                h
            };
            let mut t0 = 2000.0;
            let mut t1 = 4000.0;
            let mut h0 = hh(t0) - h_react;
            let mut h1 = hh(t1) - h_react;
            for _ in 0..60 {
                if (h1 - h0).abs() < 1e-9 {
                    break;
                }
                let t2 = t1 - h1 * (t1 - t0) / (h1 - h0);
                if !(500.0..6000.0).contains(&t2) {
                    panic!("flame secant escaped at o/f case {file}: t {t2}");
                }
                let h2 = hh(t2) - h_react;
                t0 = t1;
                h0 = h1;
                t1 = t2;
                h1 = h2;
                if h1.abs() < 1.0 {
                    break;
                }
            }
            assert!(
                (t1 - t_cea).abs() < 15.0,
                "flame temp at {file}: nufor {t1:.2} k vs cea {t_cea} k"
            );
            temps.push((o2_mol, t1));
        }
        // the sweep trend: peak near o/f 8 and falling on both sides.
        let t_stoich = temps
            .iter()
            .find(|(o, _)| (*o - 0.5).abs() < 1e-9)
            .unwrap()
            .1;
        let t_lean8 = temps
            .iter()
            .find(|(o, _)| (*o - 0.504).abs() < 1e-9)
            .unwrap()
            .1;
        let t_lean10 = temps
            .iter()
            .find(|(o, _)| (*o - 0.63).abs() < 1e-9)
            .unwrap()
            .1;
        let t_lean12 = temps
            .iter()
            .find(|(o, _)| (*o - 0.756).abs() < 1e-9)
            .unwrap()
            .1;
        let t_rich7 = temps
            .iter()
            .find(|(o, _)| (*o - 0.441).abs() < 1e-9)
            .unwrap()
            .1;
        let t_rich55 = temps
            .iter()
            .find(|(o, _)| (*o - 0.3465).abs() < 1e-9)
            .unwrap()
            .1;
        let t_rich4 = temps
            .iter()
            .find(|(o, _)| (*o - 0.252).abs() < 1e-9)
            .unwrap()
            .1;
        assert!(t_lean10 < t_lean8, "lean side must fall past o/f 8");
        assert!(t_lean12 < t_lean10, "lean side must keep falling");
        assert!(t_rich7 < t_stoich, "rich side must fall from the peak");
        assert!(t_rich55 < t_rich7, "rich side must keep falling");
        assert!(t_rich4 < t_rich55, "rich side must fall to o/f 4");
    }

    /// mole to mass fraction helper for the gate (gate-local; the modules
    /// proper work in mass fractions everywhere else).
    fn mass_from_x(x: &[f64; NS]) -> [f64; NS] {
        let mmix: f64 = (0..NS).map(|s| x[s] * M[s]).sum();
        core::array::from_fn(|s| x[s] * M[s] / mmix)
    }
}

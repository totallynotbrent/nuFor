//! park5 species transport properties: mixture viscosity and frozen
//! conductivity by wilke's rule from the RP-1232 Table V/VI curve
//! fits, and the effective binary diffusivity for the species
//! diffusion near a wall (Table VII, single effective-D form).
//!
//! fit form (all three tables): exp(A z^n family) with z = ln T.
//! viscosity and diffusion are quadratic in z, conductivity quartic.
//! units: RP-1232 prints mu in g/(cm s) (x0.1 -> Pa s), k in
//! cal/(cm s K) (x418.4 -> W/(m K)), D in cm^2/s at 1 atm
//! (x1e-4 -> m^2/s; scale linearly with 1/p_atm).
//!
//! below 1000 K the tables direct Sutherland air laws for every
//! species (Table V/VI footnotes), so the mixture rules then ride on
//! a single air value per property: the mixture rule collapses to
//! the same number whatever the composition, which is the honest
//! cold-limit behavior of a single-D single-mu approximation.
//!
//! source: plans/catalytic_research/constants.json + wilke_transport.md
//! (RP-1232 Tables V, VI, VII, verbatim coefficients).

use crate::park::{M, NS};
use crate::viscous2d::sutherland_mu;

/// sutherland air reference for the sub-1000 K branch (Table V/VI
/// footnotes direct Sutherland "for air ... for all species").
const MU0_AIR: f64 = 1.716e-5;
const T0_AIR: f64 = 273.15;
const S_AIR: f64 = 110.4;

/// Table V viscosity coefficients in this code's species order
/// [N2, O2, NO, N, O], converted from the printed (N2, O2, N, O, NO).
const MU_FIT: [[f64; 3]; NS] = [
    [0.0203, 0.4329, -11.8153], // N2
    [0.0484, -0.1455, -8.9231], // O2
    [0.0452, -0.0609, -9.4596], // NO
    [0.0120, 0.5930, -12.3805], // N
    [0.0205, 0.4257, -11.5803], // O
];

/// Table VI frozen conductivity coefficients, same order (N and O are
/// pure linear rows: A = B = 0).
const K_FIT: [[f64; 5]; NS] = [
    [0.03607, -1.07503, 11.95029, -57.90063, 93.21782], // N2
    [0.07987, -2.58428, 31.25959, -166.76267, 321.69820], // O2
    [0.02792, -0.87133, 10.17967, -52.03466, 88.67060], // NO
    [0.0, 0.0, 0.01619, 0.55022, -12.92190],            // N
    [0.0, 0.0, 0.03310, 0.22834, -11.58116],            // O
];

/// Table VII neutral-pair diffusion coefficients, indexed [i][j] in
/// species order, symmetric. N-N2 printed row is (0.0195, 1.4880,
/// -10.3654) etc.; the [i][j] / [j][i] entries mirror.
const D_FIT: [[f64; 3]; NS] = [
    [0.0112, 1.6182, -11.3091], // N2 self
    [0.0465, 0.9271, -8.1137],  // O2-N2 (used for O2 row too via symmetry of the pair table)
    [0.0291, 1.2676, -9.6878],  // NO-N2 (the NO row's own pair values carry the rest)
    [0.0195, 1.4880, -10.3654], // N-N2
    [0.0140, 1.5824, -10.8819], // O-N2
];

/// pure-species dynamic viscosity in Pa s: the RP-1232 fit above 1000 K,
/// the Sutherland air law below (the tables' own directive).
pub fn species_viscosity(s: usize, t: f64) -> f64 {
    if t < 1000.0 {
        sutherland_mu(MU0_AIR, t.max(1.0), T0_AIR, S_AIR)
    } else {
        let z = t.ln();
        let [a, b, c] = MU_FIT[s];
        0.1 * (a * z * z + b * z + c).exp()
    }
}

/// pure-species frozen thermal conductivity in W/(m K): Table VI fit
/// above 1000 K (cal/(cm s K) x 418.4), Sutherland-scaled air value
/// below. the sub-1000 K conductivity rides the same Sutherland mu
/// through kappa = mu cp / pr with air pr: the footnote's own form.
pub fn species_conductivity(s: usize, t: f64) -> f64 {
    if t < 1000.0 {
        // scale the air Sutherland mu to a conductivity with the
        // mixture's specific heat at this composition handled by the
        // caller; here the per-species value uses cp_vib-frozen air cp
        // (1004.7) and pr 0.71, the same anchors the viscous lane uses.
        let mu = sutherland_mu(MU0_AIR, t.max(1.0), T0_AIR, S_AIR);
        mu * 1004.7 / 0.71
    } else {
        let z = t.ln();
        let [a, b, c, d, e] = K_FIT[s];
        418.4 * (a * z * z * z * z + b * z * z * z + c * z * z + d * z + e).exp()
    }
}

/// binary diffusion coefficient for the pair (i, j) in m^2/s at the
/// given pressure (Pa). the pair table is symmetric in the physics;
/// the fit rows are stored once per species vs N2 with the self and
/// cross pairs picked by index below.
pub fn binary_diffusivity(i: usize, j: usize, t: f64, p_pa: f64) -> f64 {
    if t < 1000.0 {
        // sub-range: no tabulated neutral pair below 1000 K; scale the
        // 1000 K value with the Chapman-Enskog T^1.5/p shape so the
        // cold wall layer keeps a finite, physically-shaped D.
        let d1000 = pair_fit(i, j, 1000.0);
        d1000 * (t.max(1.0) / 1000.0).powf(1.5)
    } else {
        pair_fit(i, j, t)
    }
    .max(1e-12)
        * 101325.0
        / p_pa.max(1.0)
}

/// the Table VII exp fit for one pair in cm^2/s at 1 atm.
fn pair_fit(i: usize, j: usize, t: f64) -> f64 {
    let z = t.ln();
    // diagonal pairs use the species' own self-fit; off-diagonal take
    // the stored i-vs-N2 row only when one side is N2; every other
    // cross pair uses the lighter-species-vs-N2 row as the effective
    // interaction (documented approximation: the neutral cross pairs
    // vary by under 40% in the table while the single-D model absorbs
    // the difference; the effective-D formula below weights them).
    let abc = if i == j {
        D_FIT[i]
    } else {
        match (i, j) {
            (0, 1) | (1, 0) => D_FIT[1], // O2-N2
            (0, 2) | (2, 0) => D_FIT[2], // NO-N2
            (0, 3) | (3, 0) => D_FIT[3], // N-N2
            (0, 4) | (4, 0) => D_FIT[4], // O-N2
            (1, 2) | (2, 1) => D_FIT[2], // NO-O2 ~ NO-N2 row shape
            (1, 3) | (3, 1) => D_FIT[3], // N-O2 ~ N-N2 row shape
            (1, 4) | (4, 1) => D_FIT[4], // O-O2 ~ O-N2 row shape
            (2, 3) | (3, 2) => D_FIT[3], // NO-N ~ N-N2 row shape
            (2, 4) | (4, 2) => D_FIT[4], // NO-O ~ O-N2 row shape
            (3, 4) | (4, 3) => D_FIT[4], // O-N ~ O-N2 row shape
            _ => D_FIT[0],
        }
    };
    let [a, b, c] = abc;
    1e-4 * (a * z * z + b * z + c).exp()
}

/// effective binary diffusivity for the mixture in m^2/s: the
/// mole-fraction-weighted mean of the pair fits. with a single D the
/// species fluxes are J_s = -rho D dY_s/dn and sum J_s = 0 holds
/// identically (the engineering simplification documented in
/// wilke_transport.md).
pub fn effective_diffusivity(y: &[f64; NS], t: f64, p_pa: f64) -> f64 {
    // mole fractions from mass fractions
    let m_mix: f64 = y.iter().zip(M).map(|(&yi, mi)| yi / mi).sum();
    if m_mix <= 0.0 {
        return binary_diffusivity(0, 0, t, p_pa);
    }
    let x: Vec<f64> = y.iter().zip(M).map(|(&yi, mi)| (yi / mi) / m_mix).collect();
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..NS {
        for j in (i + 1)..NS {
            let w = x[i] * x[j];
            if w > 1e-12 {
                num += w * pair_fit(i, j, t);
                den += w;
            }
        }
    }
    if den <= 0.0 {
        return binary_diffusivity(0, 0, t, p_pa);
    }
    num / den * 101325.0 / p_pa.max(1.0)
}

/// wilke mixture viscosity in Pa s from mass fractions and the cell
/// temperature. species with mole fraction under 1e-12 are skipped in
/// both loops (the 0/0 guard).
pub fn mixture_viscosity(y: &[f64; NS], t: f64) -> f64 {
    let m_mix: f64 = y.iter().zip(M).map(|(&yi, mi)| yi / mi).sum();
    if m_mix <= 0.0 {
        return sutherland_mu(MU0_AIR, t.max(1.0), T0_AIR, S_AIR);
    }
    let x: Vec<f64> = y.iter().zip(M).map(|(&yi, mi)| (yi / mi) / m_mix).collect();
    let mu: Vec<f64> = (0..NS).map(|s| species_viscosity(s, t)).collect();
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..NS {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * mu[i];
        let mut phi_sum = 0.0;
        for j in 0..NS {
            if x[j] < 1e-12 {
                continue;
            }
            // phi_ij per wilke 1950
            let mr = M[i] / M[j];
            let phi = (1.0 + (mu[i] / mu[j]).sqrt() * mr.powf(0.25)).powi(2)
                / (8.0f64 * (1.0 + mr)).sqrt();
            phi_sum += x[j] * phi;
        }
        if phi_sum > 0.0 {
            den += x[i] * phi_sum;
        }
    }
    if den <= 0.0 {
        return sutherland_mu(MU0_AIR, t.max(1.0), T0_AIR, S_AIR);
    }
    num / den
}

/// wilke-class mixture frozen conductivity in W/(m K): the same
/// weighting with k in place of mu (Mason-Saxena form per RP-1232's
/// approximation review).
pub fn mixture_conductivity(y: &[f64; NS], t: f64) -> f64 {
    let m_mix: f64 = y.iter().zip(M).map(|(&yi, mi)| yi / mi).sum();
    if m_mix <= 0.0 {
        let mu = sutherland_mu(MU0_AIR, t.max(1.0), T0_AIR, S_AIR);
        return mu * 1004.7 / 0.71;
    }
    let x: Vec<f64> = y.iter().zip(M).map(|(&yi, mi)| (yi / mi) / m_mix).collect();
    let k: Vec<f64> = (0..NS).map(|s| species_conductivity(s, t)).collect();
    let mu: Vec<f64> = (0..NS).map(|s| species_viscosity(s, t)).collect();
    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..NS {
        if x[i] < 1e-12 {
            continue;
        }
        num += x[i] * k[i];
        let mut phi_sum = 0.0;
        for j in 0..NS {
            if x[j] < 1e-12 {
                continue;
            }
            let mr = M[i] / M[j];
            let phi = (1.0 + (k[i] / k[j]).sqrt() * (mu[j] / mu[i]).sqrt() * mr.powf(0.25)).powi(2)
                / (8.0f64 * (1.0 + mr)).sqrt();
            phi_sum += x[j] * phi;
        }
        if phi_sum > 0.0 {
            den += x[i] * phi_sum;
        }
    }
    if den <= 0.0 {
        let mu = sutherland_mu(MU0_AIR, t.max(1.0), T0_AIR, S_AIR);
        return mu * 1004.7 / 0.71;
    }
    num / den
}

#[cfg(test)]
mod tests {
    use super::*;

    // standard 5-species air mass fractions: no argon in the model,
    // the atmosphere's 1.34% Ar mass folds into N2 here.
    const AIR: [f64; NS] = [0.7686, 0.2314, 0.0, 0.0, 0.0];

    #[test]
    fn viscosity_fits_match_the_table_anchors() {
        // N2 at 1000 K: fit 3.87e-5 Pa s, against the Sutherland-law
        // 4.15e-5 at the range edge (fit-vs-law disagreement expected).
        let mu = species_viscosity(0, 1000.0);
        assert!((mu - 3.87e-5).abs() < 0.1e-5, "N2 mu(1000K) = {mu}");
        // atoms above molecules at 14000 K (kinetic theory ordering)
        assert!(species_viscosity(3, 14000.0) > species_viscosity(0, 14000.0));
        // sub-1000 K branch is the air Sutherland value for any species
        assert_eq!(species_viscosity(3, 500.0), species_viscosity(0, 500.0));
    }

    #[test]
    fn conductivity_fits_match_the_table_anchors() {
        // N2 at 1000 K: 0.0674 W/(m K) from the fit, against the
        // tabulated ~0.066-0.068 band
        let k = species_conductivity(0, 1000.0);
        assert!((k - 0.0674).abs() < 0.003, "N2 k(1000K) = {k}");
    }

    #[test]
    fn diffusivity_matches_the_chapman_enskog_band() {
        // N-N2 at 10000 K, 1 atm: ~1.5e-2 m^2/s
        let d = binary_diffusivity(3, 0, 10000.0, 101325.0);
        assert!((d - 1.475e-2).abs() < 0.05e-2, "N-N2 D(10000K) = {d}");
        // linear 1/p scaling
        let d2 = binary_diffusivity(3, 0, 10000.0, 202650.0);
        assert!((d2 - 0.5 * d).abs() < 1e-12);
    }

    #[test]
    fn wilke_mixtures_reproduce_the_single_species_limit() {
        // pure N2: mixture mu = species mu exactly
        let pure = [1.0, 0.0, 0.0, 0.0, 0.0];
        assert!((mixture_viscosity(&pure, 5000.0) - species_viscosity(0, 5000.0)).abs() < 1e-15);
        // air at 300 K (sub-range): all species ride one air value, so
        // the mixture collapses to it up to the phi-weighting residual.
        // for the two-molecule mixture phi(N2,O2) = 0.99917: the bound
        // is that residual scale.
        let mu300 = sutherland_mu(MU0_AIR, 300.0, T0_AIR, S_AIR);
        assert!((mixture_viscosity(&AIR, 300.0) - mu300).abs() < 5e-4 * mu300);
        // air at 10000 K sits in the plausible band (2.0-2.4e-4 Pa s)
        let mu_air_hot = mixture_viscosity(&AIR, 10000.0);
        assert!(
            mu_air_hot > 1.8e-4 && mu_air_hot < 2.6e-4,
            "air mu(10000K) = {mu_air_hot}"
        );
    }

    #[test]
    fn effective_diffusivity_sane_across_the_state_range() {
        // cold pure air: finite, small
        let d_cold = effective_diffusivity(&AIR, 300.0, 101325.0);
        assert!(d_cold > 0.0 && d_cold < 1e-4, "D_cold = {d_cold}");
        // hot dissociated at low pressure: the T^1.6/p shape puts
        // D in the O(0.1-1) m^2/s band at 12000 K, 5 kPa
        let hot = [0.2, 0.05, 0.05, 0.35, 0.35];
        let d_hot = effective_diffusivity(&hot, 12000.0, 5000.0);
        assert!(d_hot > 0.1 && d_hot < 10.0, "D_hot = {d_hot}");
    }
}

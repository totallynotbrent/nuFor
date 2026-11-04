//! equilibrium constants for GRI-Mech 3.0, derived from the NASA-7
//! thermodynamic polynomials.
//!
//! the reversible rate machinery needs the concentration-form constant
//! for every reaction:
//!
//! ```text
//! kc = exp(-dg0 / (r t)) * (p0 / (r t)) ** dn
//! ```
//!
//! `dg0` is the standard-state gibbs energy difference between products
//! and reactants at `p0 = 1 atm`, and `dn` is the change in molecule
//! count. the gibbs energy comes from the same NASA-7 fits the
//! enthalpy uses, so the constant inherits the polynomials' accuracy
//! rather than carrying a second, independent fit.
//!
//! two unit conventions meet here, and the difference matters:
//!
//! - the basis cantera reports is the classical one,
//!   `kp * (1 / (r_latm * t)) ** dn`, with concentrations in mol/l.
//! - the solver works in mol/m3, which is the same number as kmol/m3,
//!   and a concentration ten thousand times bigger than mol/l.
//!
//! so the solver's constant is cantera's times `1000 ** dn`. the
//! factor is exact, not a fudge: `p0 / (r_si * t)` in kmol/m3 and
//! `1 / (r_latm * t)` in mol/l are the same number, because
//! `r_si / p0` is `r_latm` in litres.

use crate::gri30::g_rt;

/// the gas constant in litre atmosphere per mole kelvin, the classical
/// chemical-thermodynamics value.
pub const R_LATM: f64 = 0.0820573660809596;

/// the conversion between the mole-per-litre basis and the
/// mole-per-cubic-metre basis the solver uses.
pub const MOL_PER_L_TO_MOL_PER_M3: f64 = 1000.0;

/// the change in molecule count of a reaction.
pub fn delta_nu(reactants: &[(usize, f64)], products: &[(usize, f64)]) -> f64 {
    let products: f64 = products.iter().map(|&(_, c)| c).sum();
    let reactants: f64 = reactants.iter().map(|&(_, c)| c).sum();
    products - reactants
}

/// the natural log of the equilibrium constant in the mole-per-litre
/// basis cantera reports.
///
/// the gibbs energies are carried as `g / (r t)`, so the dimensionless
/// standard-state difference is simply the weighted sum of the
/// species' `g_rt`, and the constant follows.
pub fn ln_kc_mol_per_litre(
    reactants: &[(usize, f64)],
    products: &[(usize, f64)],
    temperature: f64,
) -> f64 {
    let dn = delta_nu(reactants, products);
    let mut dg_rt = 0.0;
    for &(species, coefficient) in products {
        dg_rt += coefficient * g_rt(species, temperature);
    }
    for &(species, coefficient) in reactants {
        dg_rt -= coefficient * g_rt(species, temperature);
    }
    -dg_rt + dn * (-(R_LATM * temperature).ln())
}

/// the natural log of the equilibrium constant in the mole-per-cubic-
/// metre basis the solver uses.
///
/// this is cantera's constant shifted by `dn * ln(1000)`, which is the
/// conversion between the two concentration units raised to the
/// reaction's molecule-count change.
pub fn ln_kc_mol_per_cubic_metre(
    reactants: &[(usize, f64)],
    products: &[(usize, f64)],
    temperature: f64,
) -> f64 {
    let dn = delta_nu(reactants, products);
    ln_kc_mol_per_litre(reactants, products, temperature) + dn * MOL_PER_L_TO_MOL_PER_M3.ln()
}

/// the equilibrium constant in the mole-per-cubic-metre basis.
pub fn equilibrium_constant(
    reactants: &[(usize, f64)],
    products: &[(usize, f64)],
    temperature: f64,
) -> f64 {
    ln_kc_mol_per_cubic_metre(reactants, products, temperature).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gri30_eq_fixture;
    use crate::gri30_generated::REACTIONS;

    #[test]
    fn the_constant_matches_cantera_on_every_sampled_reaction() {
        let mut worst: f64 = 0.0;
        let mut worst_at = String::new();
        let mut compared = 0usize;
        for sample in gri30_eq_fixture::SAMPLES {
            let reaction = &REACTIONS[sample.reaction];
            let reactants: Vec<(usize, f64)> = reaction
                .reactants
                .iter()
                .map(|entry| (entry.species, entry.coefficient))
                .collect();
            let products: Vec<(usize, f64)> = reaction
                .products
                .iter()
                .map(|entry| (entry.species, entry.coefficient))
                .collect();
            let mine = ln_kc_mol_per_litre(&reactants, &products, sample.temperature);
            let error = (mine - sample.ln_kc).abs();
            compared += 1;
            if error > worst {
                worst = error;
                worst_at = format!("reaction {} at {} K", sample.reaction, sample.temperature);
            }
        }
        assert!(
            compared > 1000,
            "expected a full sweep of the mechanism, compared {}",
            compared
        );
        assert!(worst < 1.0e-9, "{} worst error {}", worst_at, worst);
    }

    #[test]
    fn the_two_unit_bases_differ_by_exactly_the_expected_shift() {
        let reaction = &REACTIONS[0];
        let reactants: Vec<(usize, f64)> = reaction
            .reactants
            .iter()
            .map(|entry| (entry.species, entry.coefficient))
            .collect();
        let products: Vec<(usize, f64)> = reaction
            .products
            .iter()
            .map(|entry| (entry.species, entry.coefficient))
            .collect();
        let dn = delta_nu(&reactants, &products);
        let litre = ln_kc_mol_per_litre(&reactants, &products, 1500.0);
        let cubic = ln_kc_mol_per_cubic_metre(&reactants, &products, 1500.0);
        let expected = dn * 1000.0f64.ln();
        assert!((cubic - litre - expected).abs() < 1.0e-12);
    }
}

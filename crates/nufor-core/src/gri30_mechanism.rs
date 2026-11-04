//! a small rust interface over the generated gri-mech 3.0 tables.
//!
//! the tables themselves live in `gri30_generated.rs`, which the
//! generator writes from cantera's own copy of the mechanism yaml;
//! this module only exposes lookups and evaluates the forward rate
//! constants.
//!
//! unit basis: concentrations are kmol/m3, which is numerically mol/l,
//! the molar gas constant is 8314.46261815324 j/(kmol k), activation
//! energies are j/kmol, and rate constants are m3/(kmol s) for
//! bimolecular laws, m6/(kmol2 s) for trimolecular ones, and 1/s for
//! unimolecular ones.
//!
//! the third-body bath follows cantera: every species contributes its
//! concentration times its collision efficiency, and a species that
//! the mechanism leaves out of the efficiency table contributes at the
//! default efficiency of one.
//!
//! every falloff reaction goes through `crate::troe`, which supplies
//! the broadening factor for the troe set and the plain lindemann
//! blend for the rest.

use crate::gri30_generated::{
    Arrhenius, Efficiency, FalloffKind, Reaction, ReactionKind, Species, Stoichiometry,
    TroeParameters, GAS_CONSTANT, REACTIONS, SOURCE_SHA256, SPECIES, SPECIES_COUNT,
};
use crate::troe::{self, TroeParams};

/// the number of species in the mechanism.
pub fn species_count() -> usize {
    SPECIES_COUNT
}

/// the number of reactions in the mechanism.
pub fn reaction_count() -> usize {
    REACTIONS.len()
}

/// the sha256 of the mechanism yaml the tables came from.
pub fn source_sha256() -> &'static str {
    SOURCE_SHA256
}

/// all species of the mechanism, in mechanism order.
pub fn species() -> &'static [Species] {
    SPECIES
}

/// the species table entry at `index`.
pub fn species_at(index: usize) -> Option<&'static Species> {
    SPECIES.get(index)
}

/// the index of the species called `name`.
pub fn species_index(name: &str) -> Option<usize> {
    SPECIES.iter().position(|s| s.name == name)
}

/// all reactions of the mechanism, in mechanism order.
pub fn reactions() -> &'static [Reaction] {
    REACTIONS
}

/// the reaction at `index`.
pub fn reaction(index: usize) -> Option<&'static Reaction> {
    REACTIONS.get(index)
}

/// the indices of the reactions that carry a falloff rate law.
pub fn falloff_reaction_indices() -> Vec<usize> {
    REACTIONS
        .iter()
        .filter(|r| matches!(r.kind, ReactionKind::Falloff))
        .map(|r| r.index)
        .collect()
}

/// the troe parameters of a falloff reaction, when it has any.
pub fn troe_parameters(index: usize) -> Option<TroeParams> {
    reaction(index)?
        .troe
        .as_ref()
        .map(|t: &TroeParameters| TroeParams {
            a: t.a,
            t3: t.t3,
            t1: t.t1,
            t2: t.t2,
        })
}

/// evaluate one arrhenius rate law at `temperature`.
pub fn arrhenius_rate(law: &Arrhenius, temperature: f64) -> f64 {
    law.a * temperature.powf(law.n) * (-law.activation_energy / (GAS_CONSTANT * temperature)).exp()
}

/// the efficiency of one species for one reaction.
///
/// a species that the mechanism leaves out of the efficiency table
/// contributes at the reaction's default efficiency.
pub fn efficiency(reaction: &Reaction, species: usize) -> f64 {
    reaction
        .efficiencies
        .iter()
        .find(|e: &&Efficiency| e.species == species)
        .map(|e| e.efficiency)
        .unwrap_or(reaction.default_efficiency)
}

/// the efficiency-weighted third-body concentration of a reaction.
///
/// `densities` holds the molar concentration of every species in
/// kmol/m3. a generic collider's default efficiency of one lets every
/// species join the bath, while a named collider's default efficiency
/// of zero leaves only the collider itself in it.
pub fn third_body_concentration(reaction: &Reaction, densities: &[f64]) -> f64 {
    let mut bath = 0.0;
    for (index, density) in densities.iter().enumerate() {
        bath += efficiency(reaction, index) * density;
    }
    bath
}

/// the reduced pressure of a falloff reaction at `temperature` and `bath`.
///
/// cantera forms it as the bath concentration times the ratio of the
/// low-pressure to the high-pressure limiting rate.
pub fn reduced_pressure(reaction: &Reaction, bath: f64, temperature: f64) -> f64 {
    let low = reaction.low.expect("falloff reaction without a low branch");
    let high = reaction
        .high
        .expect("falloff reaction without a high branch");
    let k_low = arrhenius_rate(&low, temperature);
    let k_high = arrhenius_rate(&high, temperature);
    if k_high == 0.0 {
        return 0.0;
    }
    bath * k_low / (k_high + troe::SMALL_NUMBER)
}

/// the broadening factor of a falloff reaction.
///
/// the troe set goes through `crate::troe`, and the lindemann set has
/// a factor of one.
pub fn broadening_factor(reaction: &Reaction, temperature: f64, reduced: f64) -> f64 {
    match reaction.troe {
        Some(parameters) => {
            let params = TroeParams {
                a: parameters.a,
                t3: parameters.t3,
                t1: parameters.t1,
                t2: parameters.t2,
            };
            troe::broadening_factor(&params, temperature, reduced)
        }
        None => 1.0,
    }
}

/// the forward rate constant of a reaction at `temperature`.
///
/// `densities` holds the molar concentration of every species in
/// kmol/m3. the returned constant carries the units of the rate law:
/// 1/s for a unimolecular one, m3/(kmol s) for a bimolecular one, and
/// m6/(kmol2 s) for a trimolecular one.
///
/// a falloff reaction returns its blended constant, which already
/// carries the third-body dependence through the reduced pressure.
pub fn forward_rate_constant(index: usize, temperature: f64, densities: &[f64]) -> Option<f64> {
    let reaction = reaction(index)?;
    match reaction.kind {
        ReactionKind::Arrhenius => {
            let law = reaction.forward?;
            Some(arrhenius_rate(&law, temperature))
        }
        ReactionKind::ThreeBody => {
            let law = reaction.forward?;
            let bath = third_body_concentration(reaction, densities);
            Some(arrhenius_rate(&law, temperature) * bath)
        }
        ReactionKind::Falloff => {
            let high = reaction.high?;
            let k_high = arrhenius_rate(&high, temperature);
            if k_high == 0.0 {
                return Some(0.0);
            }
            let bath = third_body_concentration(reaction, densities);
            let reduced = reduced_pressure(reaction, bath, temperature);
            let blending = reduced / (1.0 + reduced);
            let factor = match reaction.falloff {
                Some(FalloffKind::Troe) => {
                    let parameters = reaction.troe?;
                    let params = TroeParams {
                        a: parameters.a,
                        t3: parameters.t3,
                        t1: parameters.t1,
                        t2: parameters.t2,
                    };
                    troe::broadening_with_blending(Some(&params), blending, temperature, reduced)
                }
                Some(FalloffKind::Lindemann) => {
                    troe::broadening_with_blending(None, blending, temperature, reduced)
                }
                None => blending,
            };
            Some(k_high * factor)
        }
    }
}

/// the reactant list of a reaction, in the order the equation writes them.
pub fn reactants(index: usize) -> Option<&'static [Stoichiometry]> {
    reaction(index)?.reactants.get(..)
}

/// the product list of a reaction, in the order the equation writes them.
pub fn products(index: usize) -> Option<&'static [Stoichiometry]> {
    reaction(index)?.products.get(..)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::PathBuf;

    /// the oracle fixture written by the generator.
    const ORACLE_PATH: &str = "data/reference/gri30_falloff_oracle.csv";

    /// the largest relative deviation the fixture tolerates.
    const TOLERANCE: f64 = 1e-9;

    /// one parsed row of the oracle fixture.
    struct Row {
        reaction: usize,
        equation: String,
        temperature: f64,
        pressure: f64,
        mole_fractions: Vec<(usize, f64)>,
        bath: f64,
        reduced_pressure: f64,
        factor: f64,
        k_cantera: f64,
        k_projected: f64,
    }

    /// parse the fixture with plain standard-library splitting.
    fn load_rows() -> Vec<Row> {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(ORACLE_PATH);
        let text = fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("cannot read {}: {}", path.display(), e));
        let mut lines = text.lines();
        let header = lines.next().expect("fixture without a header");
        let expected = [
            "reaction_index",
            "equation",
            "temperature_K",
            "pressure_Pa",
            "mole_fraction_csv",
            "efficiency_weighted_M_mol_per_L",
            "cantera_third_body_M_mol_per_L",
            "reduced_pressure_pr",
            "broadening_factor_cantera_clean",
            "k_fwd_cantera_SI",
            "k_fwd_projected_SI",
        ];
        assert_eq!(
            header.split(',').collect::<Vec<_>>(),
            expected.to_vec(),
            "the fixture header changed"
        );
        let mut rows = Vec::new();
        for (number, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split(',').collect();
            assert_eq!(fields.len(), expected.len(), "bad row {}", number + 2);
            let mole_fractions = fields[4]
                .split(';')
                .filter(|piece| !piece.trim().is_empty())
                .map(|piece| {
                    let (name, value) = piece.split_once(':').expect("malformed mixture");
                    let index =
                        species_index(name).unwrap_or_else(|| panic!("unknown species {}", name));
                    (
                        index,
                        value.parse::<f64>().expect("malformed mole fraction"),
                    )
                })
                .collect();
            rows.push(Row {
                reaction: fields[0]
                    .parse::<usize>()
                    .expect("malformed reaction index"),
                equation: fields[1].to_string(),
                temperature: fields[2].parse::<f64>().expect("malformed temperature"),
                pressure: fields[3].parse::<f64>().expect("malformed pressure"),
                mole_fractions,
                bath: fields[5].parse::<f64>().expect("malformed bath"),
                reduced_pressure: fields[7]
                    .parse::<f64>()
                    .expect("malformed reduced pressure"),
                factor: fields[8].parse::<f64>().expect("malformed factor"),
                k_cantera: fields[9].parse::<f64>().expect("malformed cantera rate"),
                k_projected: fields[10].parse::<f64>().expect("malformed projected rate"),
            });
        }
        rows
    }

    /// the molar densities in kmol/m3 implied by a fixture row.
    fn densities(row: &Row) -> Vec<f64> {
        let total = row.pressure / (GAS_CONSTANT * row.temperature);
        let mut values = vec![0.0; species_count()];
        for (index, fraction) in &row.mole_fractions {
            values[*index] = total * fraction;
        }
        values
    }

    /// the relative difference between two positive numbers.
    fn relative(a: f64, b: f64) -> f64 {
        (a - b).abs() / b.abs().max(f64::MIN_POSITIVE)
    }

    /// every falloff reaction matches cantera's own broadening factor.
    #[test]
    fn falloff_broadening_matches_cantera() {
        let rows = load_rows();
        assert!(!rows.is_empty());
        let mut checked = 0;
        let mut worst = 0.0f64;
        let mut worst_row: Option<&Row> = None;
        for row in &rows {
            let reaction = reaction(row.reaction).expect("missing reaction");
            assert_eq!(
                reaction.equation, row.equation,
                "the equation of reaction {} moved",
                row.reaction
            );
            assert!(
                matches!(reaction.kind, ReactionKind::Falloff),
                "reaction {} is not a falloff",
                row.reaction
            );

            let densities = densities(row);
            let bath = third_body_concentration(reaction, &densities);
            let deviation = relative(bath, row.bath);
            assert!(
                deviation < 1e-12,
                "the third-body bath of reaction {} differs by {}",
                row.reaction,
                deviation
            );

            let reduced = reduced_pressure(reaction, bath, row.temperature);
            let deviation = relative(reduced, row.reduced_pressure);
            assert!(
                deviation < 1e-12,
                "the reduced pressure of reaction {} differs by {}",
                row.reaction,
                deviation
            );

            let factor = broadening_factor(reaction, row.temperature, reduced);
            let deviation = relative(factor, row.factor);
            if deviation > worst {
                worst = deviation;
                worst_row = Some(row);
            }
            assert!(
                deviation < TOLERANCE,
                "reaction {} ({}) broadens by {:.17} instead of {:.17} at {} K \
                 and a reduced pressure of {:.6}",
                row.reaction,
                row.equation,
                factor,
                row.factor,
                row.temperature,
                reduced
            );
            checked += 1;
        }
        println!(
            "checked {} broadening factors, worst relative deviation {:.3e} at {:?}",
            checked,
            worst,
            worst_row.map(|r| r.reaction)
        );
    }

    /// the blended forward rate constants match cantera's own values.
    #[test]
    fn falloff_forward_rates_match_cantera() {
        let rows = load_rows();
        let mut checked = 0;
        let mut worst = 0.0f64;
        let mut worst_description = String::new();
        for row in &rows {
            let densities = densities(row);
            let rate = forward_rate_constant(row.reaction, row.temperature, &densities)
                .expect("missing reaction");
            let deviation = relative(rate, row.k_cantera);
            if deviation > worst {
                worst = deviation;
                worst_description = format!(
                    "reaction {} ({}) at {} K",
                    row.reaction, row.equation, row.temperature
                );
            }
            assert!(
                deviation < TOLERANCE,
                "reaction {} ({}) forwards at {:.17} instead of {:.17} at {} K",
                row.reaction,
                row.equation,
                rate,
                row.k_cantera,
                row.temperature
            );
            // the fixture also carries the projected rate cantera's own
            // branches give, which must agree with its forward rate.
            let projection_deviation = relative(row.k_projected, row.k_cantera);
            assert!(
                projection_deviation < 1e-9,
                "the projected rate of reaction {} differs from cantera's by {}",
                row.reaction,
                projection_deviation
            );
            checked += 1;
        }
        println!(
            "checked {} forward rates, worst relative deviation {:.3e} at {}",
            checked, worst, worst_description
        );
    }

    /// the fixture covers every falloff reaction of the mechanism.
    #[test]
    fn fixture_covers_every_falloff_reaction() {
        let falloff = falloff_reaction_indices();
        assert_eq!(falloff.len(), 29);
        let rows = load_rows();
        let mut seen = std::collections::BTreeSet::new();
        for row in &rows {
            seen.insert(row.reaction);
        }
        for index in &falloff {
            assert!(
                seen.contains(index),
                "the fixture misses falloff reaction {}",
                index
            );
        }
        assert_eq!(seen.len(), falloff.len());

        // the conditions span the required ranges and include a
        // mixture with a non-default water efficiency and one with a
        // species outside the efficiency table.
        let mut saw_water = false;
        let mut saw_absent = false;
        for row in &rows {
            assert!(row.temperature >= 600.0 && row.temperature <= 3000.0);
            assert!(row.pressure >= 0.1 * 101_325.0 * 0.999);
            assert!(row.pressure <= 10.0 * 101_325.0 * 1.001);
            let reaction = reaction(row.reaction).expect("missing reaction");
            let water = species_index("H2O").expect("no water in the mechanism");
            if efficiency(reaction, water) != 1.0 {
                saw_water = true;
            }
            let argon = species_index("AR").expect("no argon in the mechanism");
            if efficiency(reaction, argon) == 1.0 {
                saw_absent = true;
            }
        }
        assert!(saw_water, "no case exercises a non-default efficiency");
        assert!(
            saw_absent,
            "no case exercises the default efficiency of one"
        );
    }

    /// the species and reaction bookkeeping is coherent.
    #[test]
    fn tables_are_self_consistent() {
        assert_eq!(species_count(), 53);
        assert_eq!(reaction_count(), 325);
        assert_eq!(SPECIES.len(), species_count());
        assert_eq!(REACTIONS.len(), reaction_count());
        for (index, reaction) in REACTIONS.iter().enumerate() {
            assert_eq!(reaction.index, index);
            for item in reaction.reactants.iter().chain(reaction.products.iter()) {
                assert!(item.species < species_count());
            }
        }
        assert_eq!(species_index("H2"), Some(0));
        assert_eq!(species_index("NOT_A_SPECIES"), None);
    }
}

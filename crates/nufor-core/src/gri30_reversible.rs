//! reversible rate constants and net rates of progress for GRI-Mech 3.0.
//!
//! the forward constants come from the generated mechanism tables, and
//! the backward constant follows from the equilibrium constant:
//!
//! ```text
//! kb = kf / kc
//! ```
//!
//! with `kc` the concentration-form equilibrium constant derived in
//! `gri30_equilibrium` from the NASA-7 polynomials.
//!
//! ## the concentration basis
//!
//! everything here works in kmol/m3. the solver's older hydrogen path
//! writes mol/m3, but the two are the same number: one kmol spread over
//! one cubic metre is the same amount of substance as one mol spread
//! over one litre, and kmol/m3 and mol/m3 differ by a factor of a
//! thousand, so the numerical value of a concentration is identical in
//! the two spellings only when the amount of substance is read in the
//! same unit. the generated tables carry the unit in their rate laws
//! and in the mechanisms' own comments, and this module follows them
//! rather than introducing a second convention.
//!
//! cantera's equilibrium constants use mol/l, which is numerically the
//! same as kmol/m3, so the constant transfers without a conversion
//! factor. `gri30_equilibrium` documents the conversion it does apply
//! for anyone who needs the mol/m3 spelling.
//!
//! ## falloff reactions
//!
//! a falloff reaction's forward constant already carries its pressure
//! dependence through the reduced pressure, so its backward constant is
//! obtained from the same blended forward constant. the blend commutes
//! with detailed balance: dividing the blended constant by the
//! pressure-independent equilibrium constant preserves the relation
//! `kf / kb = kc` exactly.

use crate::gri30_equilibrium;
use crate::gri30_generated::{Reaction, ReactionKind};
use crate::gri30_mechanism;

/// the equilibrium constant of a reaction, in the kmol/m3 basis.
pub fn equilibrium_constant(reaction: &Reaction, temperature: f64) -> f64 {
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
    gri30_equilibrium::ln_kc_mol_per_litre(&reactants, &products, temperature).exp()
}

/// the backward rate constant of a reaction.
///
/// the argument is the forward constant, so a falloff reaction's
/// pressure-blended value is the one that gets divided.
pub fn backward_rate_constant(reaction: &Reaction, forward: f64, temperature: f64) -> f64 {
    if !reaction.reversible {
        return 0.0;
    }
    let kc = equilibrium_constant(reaction, temperature);
    if kc <= 0.0 || !kc.is_finite() {
        return 0.0;
    }
    forward / kc
}

/// the net rate of progress of a reaction, in kmol/m3/s.
///
/// `densities` holds every species' concentration in kmol/m3. the
/// forward term raises the reactants to their powers and the backward
/// term lowers the products, which is the standard form for a
/// reversible elementary reaction.
pub fn net_rate_of_progress(index: usize, temperature: f64, densities: &[f64]) -> f64 {
    let reaction = match gri30_mechanism::reaction(index) {
        Some(reaction) => reaction,
        None => return 0.0,
    };
    let forward = match gri30_mechanism::forward_rate_constant(index, temperature, densities) {
        Some(value) => value,
        None => return 0.0,
    };
    let mut rate = forward;
    for entry in reaction.reactants {
        rate *= densities[entry.species].powf(entry.coefficient);
    }
    if !reaction.reversible {
        return rate;
    }
    let backward = backward_rate_constant(reaction, forward, temperature);
    let mut reverse = backward;
    for entry in reaction.products {
        reverse *= densities[entry.species].powf(entry.coefficient);
    }
    rate - reverse
}

/// the per-species production rates of the mechanism, in kmol/m3/s.
///
/// the result is indexed by the mechanism's species indices.
pub fn net_production_rates(temperature: f64, densities: &[f64], out: &mut [f64]) {
    for value in out.iter_mut() {
        *value = 0.0;
    }
    for index in 0..gri30_mechanism::reaction_count() {
        let reaction = match gri30_mechanism::reaction(index) {
            Some(reaction) => reaction,
            None => continue,
        };
        let rate = net_rate_of_progress(index, temperature, densities);
        if rate == 0.0 {
            continue;
        }
        for entry in reaction.reactants {
            out[entry.species] -= entry.coefficient * rate;
        }
        for entry in reaction.products {
            out[entry.species] += entry.coefficient * rate;
        }
    }
}

/// whether a reaction's forward constant carries a third body.
///
/// exposed so the caller can reason about the density argument without
/// reaching into the generated tables.
pub fn consumes_third_body(reaction: &Reaction) -> bool {
    matches!(
        reaction.kind,
        ReactionKind::ThreeBody | ReactionKind::Falloff
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_irreversible_reaction_has_no_backward_rate() {
        let irreversible = gri30_mechanism::reactions()
            .iter()
            .find(|reaction| !reaction.reversible)
            .expect("the mechanism has irreversible reactions");
        assert_eq!(backward_rate_constant(irreversible, 1.0e12, 1500.0), 0.0);
    }

    #[test]
    fn the_backward_constant_recovers_the_forward_one_through_the_equilibrium_constant() {
        for reaction in gri30_mechanism::reactions() {
            if !reaction.reversible {
                continue;
            }
            let forward = 1.0e8;
            let backward = backward_rate_constant(reaction, forward, 1500.0);
            let kc = equilibrium_constant(reaction, 1500.0);
            assert!(
                ((forward / backward) - kc).abs() <= 1.0e-9 * kc.abs(),
                "reaction {} broke the detailed-balance relation",
                reaction.equation
            );
        }
    }

    /// the slot of one species in a concentration vector.
    fn species_slot(name: &str) -> usize {
        gri30_mechanism::species_index(name).expect("the species is in the mechanism")
    }

    #[test]
    fn a_named_collider_joins_the_bath_alone() {
        // the hydrogen-plus-oxygen reaction with argon as its named
        // collider draws argon from the bath at an efficiency of one and
        // nothing from any other species.
        let reaction = gri30_mechanism::reactions()
            .iter()
            .find(|reaction| reaction.index == 36)
            .expect("the mechanism has the argon-collided oxidation");
        assert_eq!(reaction.default_efficiency, 0.0);
        assert_eq!(
            reaction.efficiencies.len(),
            1,
            "the argon collider carries one efficiency"
        );

        let mut densities = vec![0.0f64; gri30_mechanism::species_count()];
        let argon = species_slot("AR");
        densities[argon] = 2.0;

        let bath = gri30_mechanism::third_body_concentration(reaction, &densities);
        assert!((bath - 2.0).abs() <= 1.0e-12, "the bath held {:?}", bath);
    }

    #[test]
    fn a_generic_collider_lets_every_species_join_the_bath() {
        // a generic 'M' collider's default efficiency of one lets every
        // species contribute, including the collider species itself.
        let reaction = gri30_mechanism::reactions()
            .iter()
            .find(|reaction| {
                matches!(reaction.kind, ReactionKind::ThreeBody)
                    && reaction.default_efficiency == 1.0
                    && reaction.efficiencies.len() > 1
            })
            .expect("the mechanism has a generic third body");
        assert_eq!(reaction.default_efficiency, 1.0);

        let mut densities = vec![0.0f64; gri30_mechanism::species_count()];
        densities[0] = 1.0;
        densities[1] = 3.0;

        // the bath sums the table's efficiencies for the species they
        // name and the default efficiency for everyone else.
        let mut expected = 0.0;
        for (slot, density) in densities.iter().enumerate() {
            if *density == 0.0 {
                continue;
            }
            let efficiency = reaction
                .efficiencies
                .iter()
                .find(|entry| entry.species == slot)
                .map(|entry| entry.efficiency)
                .unwrap_or(reaction.default_efficiency);
            expected += efficiency * density;
        }

        let bath = gri30_mechanism::third_body_concentration(reaction, &densities);
        assert!(
            (bath - expected).abs() <= 1.0e-12,
            "the bath held {:?} against an expected {:?}",
            bath,
            expected
        );
        assert!(
            expected > 0.0,
            "the chosen reaction drew no species into the bath"
        );
    }

    #[test]
    fn the_bath_multiplies_both_branches_of_the_rate() {
        // doubling the bath doubles the magnitude of the net rate, since
        // the reverse branch picks up the same bath as the forward one.
        let reaction = gri30_mechanism::reactions()
            .iter()
            .find(|reaction| reaction.index == 36)
            .expect("the mechanism has the argon-collided oxidation");
        let mut densities = vec![0.0f64; gri30_mechanism::species_count()];
        densities[species_slot("H")] = 1.0e-3;
        densities[species_slot("O2")] = 1.0e-2;
        densities[species_slot("HO2")] = 1.0e-6;
        densities[species_slot("AR")] = 1.0;

        let single = net_rate_of_progress(reaction.index, 1800.0, &densities);
        densities[species_slot("AR")] *= 2.0;
        let doubled = net_rate_of_progress(reaction.index, 1800.0, &densities);
        assert!(single.abs() > 0.0, "the reaction produced no rate to scale");
        let ratio = doubled / single;
        assert!(
            (ratio - 2.0).abs() <= 1.0e-9,
            "doubling the bath scaled the rate by {:?}",
            ratio
        );
    }

    #[test]
    fn a_zero_density_mixture_gives_a_zero_net_rate() {
        let densities = vec![0.0f64; gri30_mechanism::species_count()];
        for index in 0..gri30_mechanism::reaction_count() {
            assert_eq!(net_rate_of_progress(index, 1500.0, &densities), 0.0);
        }
    }

    #[test]
    fn the_net_rates_match_cantera() {
        // cantera's own net rates of progress, for the same mixtures, are
        // the reference. the comparison covers every reaction of the
        // expanded list.
        let mut worst_relative = 0.0f64;
        let mut worst_description = String::new();
        for (case_index, mixture) in crate::gri30_netrate_fixture::MIXTURES.iter().enumerate() {
            let expected = &crate::gri30_netrate_fixture::NET_RATES[case_index];
            for (reaction_index, &reference) in expected.iter().enumerate() {
                let mine =
                    net_rate_of_progress(reaction_index, mixture.temperature, &mixture.densities);
                let scale = reference.abs().max(mine.abs());
                if scale == 0.0 {
                    continue;
                }
                let relative = ((mine - reference).abs() / scale).abs();
                if relative > worst_relative {
                    worst_relative = relative;
                    let reaction = gri30_mechanism::reaction(reaction_index)
                        .map(|reaction| reaction.equation.to_string())
                        .unwrap_or_else(|| "?".to_string());
                    worst_description = format!(
                        "reaction {} ({}) at {:.1} K: mine {} against cantera {}",
                        reaction_index, reaction, mixture.temperature, mine, reference
                    );
                }
            }
        }
        assert!(
            worst_relative < 1.0e-9,
            "the net rates drifted from cantera by {}: {}",
            worst_relative,
            worst_description
        );
    }

    #[test]
    fn the_production_rates_cancel_for_a_closed_system() {
        // every reaction conserves the atoms it carries, so the totals
        // of the elemental composition have to cancel exactly.
        let mut densities = vec![0.0f64; gri30_mechanism::species_count()];
        let mut rng_seed = 12345_u64;
        for value in densities.iter_mut() {
            rng_seed = rng_seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            *value = ((rng_seed >> 33) as f64) / 1000.0;
        }
        let mut rates = vec![0.0f64; gri30_mechanism::species_count()];
        net_production_rates(1500.0, &densities, &mut rates);

        let names = ["H", "C", "N", "O", "Ar"];
        for (column, element) in names.iter().enumerate() {
            let mut total = 0.0;
            let mut magnitude = 0.0;
            for (species, rate) in rates.iter().enumerate() {
                let contribution = crate::gri30_thermo::COMPOSITION[species][column] * rate;
                total += contribution;
                magnitude += contribution.abs();
            }
            // the individual contributions dwarf the sum, so the check has
            // to be relative to the largest of them.
            let relative = if magnitude > 0.0 {
                total.abs() / magnitude
            } else {
                0.0
            };
            assert!(
                relative < 1.0e-12,
                "element {} was not conserved, residual {} against a magnitude of {}",
                element,
                total,
                magnitude
            );
        }
    }
}

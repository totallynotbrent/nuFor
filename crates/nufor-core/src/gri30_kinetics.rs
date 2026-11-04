//! an allocation-free path through the gri-mech kinetics.
//!
//! the reference path in `gri30_reversible` builds two vectors per
//! reaction per call, which suits a verification sweep but not a march.
//! a zero-dimensional reactor evaluates the mechanism tens of thousands
//! of times, so the cost of the allocations dominates the cost of the
//! arithmetic.
//!
//! this module computes the same quantities from the same tables with
//! no heap traffic:
//!
//! - the third-body efficiencies resolve through a dense row per
//!   reaction, built once,
//! - the equilibrium constant accumulates straight from the static
//!   stoichiometry, without the intermediate vectors,
//! - the stoichiometric powers multiply directly for the integer
//!   coefficients the mechanism uses.
//!
//! the arithmetic follows the reference path step for step, so the two
//! agree to round-off. the module's tests hold them to each other.

use crate::gri30;
use crate::gri30_generated::{Arrhenius, Reaction, ReactionKind, Stoichiometry};
use crate::gri30_mechanism;
use crate::troe;

/// the gas constant in litre atmosphere per mole kelvin.
///
/// the equilibrium constant uses the classical standard state, which is
/// spelled in litre atmospheres.
const R_LATM: f64 = 0.0820573660809596;

/// the number of species the mechanism carries.
pub const NS: usize = crate::gri30_thermo::NS_GRI30;

/// the number of reactions the mechanism carries.
pub const NR: usize = crate::gri30_generated::REACTION_COUNT;

/// the dense efficiency rows, one per reaction.
///
/// each row holds the collisional efficiency of every species, so the
/// bath of a third-body reaction is a dot product rather than a search.
struct EfficiencyTable {
    rows: Vec<Vec<f64>>,
}

impl EfficiencyTable {
    /// build the dense rows from the static tables.
    fn build() -> Self {
        let mut rows = Vec::with_capacity(NR);
        for reaction in gri30_mechanism::reactions() {
            let mut row = vec![reaction.default_efficiency; NS];
            for entry in reaction.efficiencies {
                row[entry.species] = entry.efficiency;
            }
            rows.push(row);
        }
        EfficiencyTable { rows }
    }

    /// the bath of a reaction, the efficiency-weighted concentration.
    fn bath(&self, reaction_index: usize, concentrations: &[f64]) -> f64 {
        let row = &self.rows[reaction_index];
        let mut bath = 0.0f64;
        for (index, value) in concentrations.iter().enumerate() {
            bath += row[index] * value;
        }
        bath
    }
}

fn efficiency_table() -> &'static EfficiencyTable {
    use std::sync::OnceLock;
    static TABLE: OnceLock<EfficiencyTable> = OnceLock::new();
    TABLE.get_or_init(EfficiencyTable::build)
}

/// the natural log of a reaction's equilibrium constant.
///
/// the accumulation follows the reference path's order, so the two
/// produce the same bits.
fn ln_equilibrium_constant(reaction: &Reaction, temperature: f64) -> f64 {
    let mut reactant_moles = 0.0f64;
    let mut product_moles = 0.0f64;
    let mut dg_rt = 0.0f64;
    for entry in reaction.products {
        product_moles += entry.coefficient;
        dg_rt += entry.coefficient * gri30::g_rt(entry.species, temperature);
    }
    for entry in reaction.reactants {
        reactant_moles += entry.coefficient;
        dg_rt -= entry.coefficient * gri30::g_rt(entry.species, temperature);
    }
    let delta = product_moles - reactant_moles;
    -dg_rt + delta * (-(R_LATM * temperature).ln())
}

/// the arrhenius rate law, in the mechanism's own units.
#[inline]
fn arrhenius_rate(law: &Arrhenius, temperature: f64) -> f64 {
    law.a
        * temperature.powf(law.n)
        * (-law.activation_energy / (crate::gri30_generated::GAS_CONSTANT * temperature)).exp()
}

/// the forward rate constant of a reaction, without the stoichiometry.
fn forward_rate_constant(
    reaction: &Reaction,
    index: usize,
    temperature: f64,
    concentrations: &[f64],
) -> f64 {
    let table = efficiency_table();
    match reaction.kind {
        ReactionKind::Arrhenius => match reaction.forward {
            Some(law) => arrhenius_rate(&law, temperature),
            None => 0.0,
        },
        ReactionKind::ThreeBody => match reaction.forward {
            Some(law) => {
                let base = arrhenius_rate(&law, temperature);
                base * table.bath(index, concentrations)
            }
            None => 0.0,
        },
        ReactionKind::Falloff => {
            let high = match reaction.high {
                Some(law) => law,
                None => return 0.0,
            };
            let k_high = arrhenius_rate(&high, temperature);
            if k_high == 0.0 {
                return 0.0;
            }
            let low = match reaction.low {
                Some(law) => law,
                None => return 0.0,
            };
            let k_low = arrhenius_rate(&low, temperature);
            let bath = table.bath(index, concentrations);
            let reduced = bath * k_low / (k_high + troe::SMALL_NUMBER);
            let blending = reduced / (1.0 + reduced);
            let broadening = match reaction.troe {
                Some(parameters) => {
                    let params = troe::TroeParams {
                        a: parameters.a,
                        t3: parameters.t3,
                        t1: parameters.t1,
                        t2: parameters.t2,
                    };
                    troe::broadening_factor(&params, temperature, reduced)
                }
                None => 1.0,
            };
            k_high * broadening * blending
        }
    }
}

/// the concentration product of a stoichiometric set.
///
/// the coefficients in the mechanism are almost all one or two, and
/// multiplying directly is both faster and exact.
#[inline]
fn concentration_product(concentrations: &[f64], entries: &[Stoichiometry]) -> f64 {
    let mut product = 1.0f64;
    for entry in entries {
        let value = concentrations[entry.species];
        if entry.coefficient == 1.0 {
            product *= value;
        } else if entry.coefficient == 2.0 {
            product *= value * value;
        } else {
            product *= value.powf(entry.coefficient);
        }
    }
    product
}

/// the evaluated mechanism at one state.
#[derive(Clone)]
pub struct Evaluation {
    /// the net production rate of every species, kmol/m3/s.
    pub omega: Vec<f64>,
    /// the forward rate constant of every reaction.
    pub forward_constant: Vec<f64>,
    /// the backward rate constant of every reaction.
    pub backward_constant: Vec<f64>,
}

/// evaluate the mechanism at one temperature and composition.
///
/// the result carries the per-reaction constants the jacobian needs, so
/// the derivative assembly does not repeat the rate work.
pub fn evaluate(temperature: f64, concentrations: &[f64]) -> Evaluation {
    let mut omega = vec![0.0f64; NS];
    let mut forward_constant = vec![0.0f64; NR];
    let mut backward_constant = vec![0.0f64; NR];
    for (index, reaction) in gri30_mechanism::reactions().iter().enumerate() {
        let kf = forward_rate_constant(reaction, index, temperature, concentrations);
        forward_constant[index] = kf;
        if kf == 0.0 {
            continue;
        }
        let mut kb = 0.0f64;
        if reaction.reversible {
            let ln_kc = ln_equilibrium_constant(reaction, temperature);
            kb = kf / ln_kc.exp();
        }
        backward_constant[index] = kb;
        let forward = kf * concentration_product(concentrations, reaction.reactants);
        let backward = kb * concentration_product(concentrations, reaction.products);
        let net = forward - backward;
        if net == 0.0 {
            continue;
        }
        for entry in reaction.reactants {
            omega[entry.species] -= entry.coefficient * net;
        }
        for entry in reaction.products {
            omega[entry.species] += entry.coefficient * net;
        }
    }
    Evaluation {
        omega,
        forward_constant,
        backward_constant,
    }
}

/// the net production rates of the mechanism, kmol/m3/s.
///
/// the concentrations are molar densities in kmol/m3.
pub fn net_production_rates(temperature: f64, concentrations: &[f64], out: &mut [f64]) {
    for value in out.iter_mut() {
        *value = 0.0;
    }
    for (index, reaction) in gri30_mechanism::reactions().iter().enumerate() {
        let kf = forward_rate_constant(reaction, index, temperature, concentrations);
        if kf == 0.0 {
            continue;
        }
        let mut kb = 0.0f64;
        if reaction.reversible {
            let ln_kc = ln_equilibrium_constant(reaction, temperature);
            kb = kf / ln_kc.exp();
        }
        let forward = kf * concentration_product(concentrations, reaction.reactants);
        let backward = kb * concentration_product(concentrations, reaction.products);
        let net = forward - backward;
        if net == 0.0 {
            continue;
        }
        for entry in reaction.reactants {
            out[entry.species] -= entry.coefficient * net;
        }
        for entry in reaction.products {
            out[entry.species] += entry.coefficient * net;
        }
    }
}

/// the net rate of progress of one reaction, kmol/m3/s.
pub fn net_rate_of_progress(index: usize, temperature: f64, concentrations: &[f64]) -> f64 {
    let reaction = match gri30_mechanism::reaction(index) {
        Some(reaction) => reaction,
        None => return 0.0,
    };
    let kf = forward_rate_constant(reaction, index, temperature, concentrations);
    if kf == 0.0 {
        return 0.0;
    }
    let forward = kf * concentration_product(concentrations, reaction.reactants);
    if !reaction.reversible {
        return forward;
    }
    let ln_kc = ln_equilibrium_constant(reaction, temperature);
    let kb = kf / ln_kc.exp();
    forward - kb * concentration_product(concentrations, reaction.products)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gri30_reversible;

    /// a small deterministic generator, so the test never needs a
    /// random-number dependency.
    fn pseudo_random(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((*seed >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// the fast path reproduces the reference path's production rates.
    #[test]
    fn fast_kinetics_matches_the_reference() {
        let mut seed = 0x5eed_1234_u64;
        let mut worst: f64 = 0.0;
        for _ in 0..20 {
            let temperature = 400.0 + 2600.0 * pseudo_random(&mut seed);
            let mut concentrations = vec![0.0f64; NS];
            let mut total = 0.0f64;
            for value in concentrations.iter_mut() {
                *value = pseudo_random(&mut seed).powf(8.0) * 0.1;
                total += *value;
            }
            for value in concentrations.iter_mut() {
                *value /= total;
            }
            let mut fast = vec![0.0f64; NS];
            let mut reference = vec![0.0f64; NS];
            net_production_rates(temperature, &concentrations, &mut fast);
            gri30_reversible::net_production_rates(temperature, &concentrations, &mut reference);
            for s in 0..NS {
                let scale = fast[s].abs().max(reference[s].abs()).max(1e-30);
                let error = (fast[s] - reference[s]).abs() / scale;
                if error > worst {
                    worst = error;
                }
            }
        }
        assert!(worst < 1e-12, "worst relative difference {worst}");
    }

    /// the fast path reproduces the reference path's reaction rates.
    #[test]
    fn fast_reaction_rates_match_the_reference() {
        let mut seed = 0xabcd_9876_u64;
        let mut worst: f64 = 0.0;
        for _ in 0..10 {
            let temperature = 500.0 + 2500.0 * pseudo_random(&mut seed);
            let mut concentrations = vec![0.0f64; NS];
            let mut total = 0.0f64;
            for value in concentrations.iter_mut() {
                *value = pseudo_random(&mut seed).powf(6.0) * 0.05;
                total += *value;
            }
            for value in concentrations.iter_mut() {
                *value /= total;
            }
            for index in 0..NR {
                let fast = net_rate_of_progress(index, temperature, &concentrations);
                let reference =
                    gri30_reversible::net_rate_of_progress(index, temperature, &concentrations);
                let scale = fast.abs().max(reference.abs()).max(1e-30);
                let error = (fast - reference).abs() / scale;
                if error > worst {
                    worst = error;
                }
            }
        }
        assert!(worst < 1e-12, "worst relative difference {worst}");
    }
}

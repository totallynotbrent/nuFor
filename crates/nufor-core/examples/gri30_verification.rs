//! dumps the mechanism's verification data for the documentation figures.
//!
//! the output is a comma-separated table with one row per evaluated
//! mixture and reaction:
//!
//!     mixture,temperature_K,pressure_Pa,reaction,ours,cantera
//!
//! `ours` is the net rate of progress this crate computes and `cantera`
//! is the value in the checked-in fixture, which came from cantera
//! 3.2.0 itself. the ratio column lets the figures show the agreement
//! without re-reading rust.
//!
//! usage:
//!
//!     cargo run -p nufor-core --example gri30_verification > out.csv

use nufor_core::gri30_mechanism;
use nufor_core::gri30_netrate_fixture;
use nufor_core::gri30_reversible;

fn main() {
    println!("mixture,temperature_K,pressure_Pa,reaction,ours,cantera");
    for (mix_index, mixture) in gri30_netrate_fixture::MIXTURES.iter().enumerate() {
        for (reaction_index, cantera) in gri30_netrate_fixture::NET_RATES[mix_index]
            .iter()
            .enumerate()
        {
            let ours = gri30_reversible::net_rate_of_progress(
                reaction_index,
                mixture.temperature,
                &mixture.densities,
            );
            println!(
                "{mix_index},{:.6},{:.3},{reaction_index},{ours:.12e},{cantera:.12e}",
                mixture.temperature, mixture.pressure
            );
        }
    }

    println!();
    println!("# falloff reactions and their broadening width");
    println!("reaction,equation,a,T3,T1,T2");
    for index in gri30_mechanism::falloff_reaction_indices() {
        let reaction = &gri30_mechanism::reactions()[index];
        if let Some(parameters) = &reaction.troe {
            println!(
                "{index},{},{:.10e},{:.10e},{:.10e},{:.10e}",
                reaction.equation.replace(',', "|"),
                parameters.a,
                parameters.t3,
                parameters.t1,
                parameters.t2
            );
        }
    }

    println!();
    println!("# broadening factors of the falloff reactions");
    println!("reaction,kind,temperature_K,reduced_pressure,broadening_factor");
    let grid: [f64; 9] = [-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0];
    for index in gri30_mechanism::falloff_reaction_indices() {
        let reaction = &gri30_mechanism::reactions()[index];
        let kind = if reaction.troe.is_some() {
            "troe"
        } else {
            "plain"
        };
        for temperature in [800.0_f64, 1500.0, 2500.0] {
            for exponent in grid {
                let reduced = 10.0_f64.powf(exponent);
                let factor = gri30_mechanism::broadening_factor(reaction, temperature, reduced);
                println!("{index},{kind},{temperature:.1},{reduced:.12e},{factor:.12e}");
            }
        }
    }

    println!();
    println!("# third-body bath contribution per species");
    println!("reaction,kind,species,density_kmol_m3,efficiency,contribution");
    let mixture = &gri30_netrate_fixture::MIXTURES[2];
    for reaction in gri30_mechanism::reactions() {
        if !reaction.third_body {
            continue;
        }
        let kind = if reaction.default_efficiency == 0.0 {
            "named"
        } else {
            "generic"
        };
        for species in 0..gri30_mechanism::species_count() {
            let density = mixture.densities[species];
            let efficiency = gri30_mechanism::efficiency(reaction, species);
            if density == 0.0 && efficiency == 0.0 {
                continue;
            }
            println!(
                "{},{kind},{species},{density:.12e},{efficiency:.6},{:.12e}",
                reaction.index,
                density * efficiency
            );
        }
    }

    println!();
    println!("# third-body reactions and their default efficiency");
    println!("reaction,equation,kind,default_efficiency,efficiency_entries");
    for reaction in gri30_mechanism::reactions() {
        if !reaction.third_body {
            continue;
        }
        let kind = if reaction.default_efficiency == 0.0 {
            "named"
        } else {
            "generic"
        };
        println!(
            "{},{},{kind},{:.1},{}",
            reaction.index,
            reaction.equation.replace(',', "|"),
            reaction.default_efficiency,
            reaction.efficiencies.len()
        );
    }
}

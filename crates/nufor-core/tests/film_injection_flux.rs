//! behavioural gate for the blown-in film's mass, species and enthalpy.
//!
//! the engine wall lane's film blows THROUGH the wall into the gas, so
//! it injects mass, carries species, and arrives with its own
//! enthalpy. the assertions below are physical invariants of that
//! contract, taken from the public api only.

use nufor_core::combustion_march::{h2o2_ic_energy, CombustionState};
use nufor_core::combustion_thermo as ct;
use nufor_core::combustion_viscous_axi::{
    add_axi_combustion_viscous, axi_cell_volume, axi_combustion_transport, wall_diagnostics,
    AxiWallOpts,
};
use nufor_core::film_cooling::{coolant_composition, coolant_enthalpy_at};
use nufor_core::{grid2d, ConservedState2d, Grid2d};

/// the composition kinds the case file's
/// `physics.wall_film_injection_composition` accepts.
const CONFIG_COMPOSITION_CHOICES: [&str; 5] = ["h2", "o2", "n2", "water", "stoich"];

fn chamber() -> Grid2d {
    grid2d(8, 8, 0.0, 0.08, 0.0, 0.08).expect("the chamber's grid")
}

fn resting_gas(g: &Grid2d) -> (ConservedState2d, CombustionState) {
    let n = g.nx * g.ny;
    let rho = 12.64f64;
    let e_int = h2o2_ic_energy(rho, 7.0e6);
    (
        ConservedState2d {
            rho: vec![rho; n],
            mx: vec![0.0; n],
            my: vec![0.0; n],
            e: vec![rho * e_int; n],
        },
        CombustionState::stoich_h2o2(n),
    )
}

/// the density of the mixture's molar amount, in kmol per kilogram.
fn kmol_per_kg(y: &[f64; ct::NS]) -> Vec<f64> {
    (0..ct::NS).map(|s| y[s] / ct::M[s]).collect()
}

#[test]
fn every_configurable_film_mixture_has_mass_fractions_summing_to_one() {
    for kind in CONFIG_COMPOSITION_CHOICES {
        let y = coolant_composition(kind);
        let sum: f64 = y.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1.0e-12,
            "{kind}: mass fractions sum to {sum}, not one"
        );
        assert!(
            y.iter().all(|&f| f >= 0.0),
            "{kind}: a mass fraction went negative"
        );
        let carrying = y.iter().filter(|&&f| f > 0.0).count();
        assert!(carrying >= 1, "{kind}: the mixture carries no species");
    }
}

#[test]
fn a_hotter_film_carries_more_enthalpy_than_a_colder_one_of_the_same_composition() {
    for kind in CONFIG_COMPOSITION_CHOICES {
        let y = coolant_composition(kind);
        let cold = coolant_enthalpy_at(300.0, &y);
        let warm = coolant_enthalpy_at(600.0, &y);
        assert!(
            warm > cold,
            "{kind}: enthalpy fell from {cold} to {warm} as the film heated"
        );
        let cp = nufor_core::film_cooling::coolant_capacity_at(450.0, &y);
        let slope = (warm - cold) / 300.0;
        assert!(
            (slope - cp).abs() / cp < 0.2,
            "{kind}: the enthalpy slope {slope} disagrees with the capacity {cp}"
        );
    }
}

#[test]
fn different_compositions_at_the_same_temperature_carry_different_enthalpy() {
    let temperature = 500.0;
    let h2 = coolant_enthalpy_at(temperature, &coolant_composition("h2"));
    let o2 = coolant_enthalpy_at(temperature, &coolant_composition("o2"));
    let n2 = coolant_enthalpy_at(temperature, &coolant_composition("n2"));
    let water = coolant_enthalpy_at(temperature, &coolant_composition("water"));
    let stoich = coolant_enthalpy_at(temperature, &coolant_composition("stoich"));

    let mut seen = [
        ("h2", h2),
        ("o2", o2),
        ("n2", n2),
        ("water", water),
        ("stoich", stoich),
    ];
    seen.sort_by(|a, b| a.1.total_cmp(&b.1));
    for pair in seen.windows(2) {
        assert!(
            (pair[0].1 - pair[1].1).abs() > 1.0e3,
            "{} and {} carry the same enthalpy at {temperature} K",
            pair[0].0,
            pair[1].0
        );
    }
}

#[test]
fn the_configuration_choice_maps_to_the_expected_species_mixture() {
    let h2 = coolant_composition("h2");
    assert_eq!(h2[ct::IH2], 1.0, "\"h2\" is not pure hydrogen");
    assert!(
        h2.iter()
            .enumerate()
            .all(|(s, &f)| s == ct::IH2 || f == 0.0),
        "\"h2\" drags other species along"
    );

    let stoich = coolant_composition("stoich");
    assert!(
        stoich
            .iter()
            .enumerate()
            .all(|(s, &f)| s == ct::IH2 || s == ct::IO2 || f == 0.0),
        "\"stoich\" carries species beyond hydrogen and oxygen"
    );
    assert!(
        (stoich[ct::IH2] - 2.016 / 18.016).abs() < 1.0e-12,
        "\"stoich\" hydrogen fraction is {}",
        stoich[ct::IH2]
    );
    assert!(
        (stoich[ct::IO2] - 16.0 / 18.016).abs() < 1.0e-12,
        "\"stoich\" oxygen fraction is {}",
        stoich[ct::IO2]
    );

    // the unburned charge is stoichiometric in MOLES, not in kilograms.
    let kmol = kmol_per_kg(&stoich);
    // the charge names rounded molar masses, so the ratio lands on two
    // to five decimal places rather than exactly.
    let ratio = kmol[ct::IH2] / kmol[ct::IO2];
    assert!(
        (ratio - 2.0).abs() < 1.0e-4,
        "\"stoich\" hydrogen to oxygen molar ratio is {ratio}"
    );

    let unknown = coolant_composition("no-such-kind");
    assert_eq!(
        unknown[ct::IH2],
        1.0,
        "an unnamed coolant does not fall back to hydrogen"
    );
}

#[test]
fn the_blown_in_film_enters_the_chamber_as_mass_species_and_enthalpy() {
    let g = chamber();
    let nx = g.nx;
    let n = g.nx * g.ny;
    let (state, species) = resting_gas(&g);
    let dist = |_x: f64, y: f64| y - 0.02;
    let dt = 1.0e-3;
    let mass_flux = 5.0f64;
    let injection_temperature = 350.0f64;
    let coolant = coolant_composition("n2");

    let opts = AxiWallOpts {
        temperature: 0.0,
        emissivity: 0.0,
        film_coefficient: 0.0,
        film_injection_mass_flux: mass_flux,
        film_injection_temperature: injection_temperature,
        film_injection_composition: coolant,
        ..Default::default()
    };

    let mut advanced = state.clone();
    let mut advanced_species = species.clone();
    let transport = axi_combustion_transport(&state, &species, &g).expect("the transport snapshot");
    add_axi_combustion_viscous(
        &mut advanced,
        &mut advanced_species,
        &g,
        &transport,
        dt,
        opts,
        Some(&dist),
    )
    .expect("the wall application");

    let report =
        wall_diagnostics(&state, &species, &g, 4, Some(&dist), opts).expect("the wall report");
    assert!(
        report.film_injection_mass_flow_kg_per_second > 0.0,
        "the blown-in film crossed no area at all"
    );
    assert!(
        report.film_injection_area_square_metres > 0.0,
        "the blown-in film met no wall"
    );

    // the film's mass reached the chamber, and the report agrees.
    let mass_gain: f64 = (0..n)
        .map(|k| (advanced.rho[k] - state.rho[k]) * axi_cell_volume(&g, k % nx, k / nx))
        .sum();
    let booked = report.film_injection_mass_flow_kg_per_second * dt;
    assert!(
        mass_gain > 0.0,
        "the chamber took no mass from the blown-in film"
    );
    assert!(
        (mass_gain - booked).abs() / booked < 1.0e-9,
        "the chamber took {mass_gain} kg of coolant, the film booked {booked} kg"
    );

    // the film's enthalpy arrived with it: the energy per kilogram is
    // the mixture's enthalpy at the injection temperature.
    let energy_gain: f64 = (0..n)
        .map(|k| (advanced.e[k] - state.e[k]) * axi_cell_volume(&g, k % nx, k / nx))
        .sum();
    let expected_enthalpy = coolant_enthalpy_at(injection_temperature, &coolant);
    let measured_enthalpy = energy_gain / mass_gain;
    assert!(
        (measured_enthalpy - expected_enthalpy).abs() / expected_enthalpy.abs() < 1.0e-6,
        "the film arrived with {measured_enthalpy} J/kg, the mixture carries \
         {expected_enthalpy} J/kg at {injection_temperature} K"
    );

    // the film's species landed in the cells that took its mass.
    let mut receiver = 0usize;
    let mut receiver_gain = 0.0f64;
    for k in 0..n {
        let gain = (advanced.rho[k] - state.rho[k]) * axi_cell_volume(&g, k % nx, k / nx);
        if gain > receiver_gain {
            receiver_gain = gain;
            receiver = k;
        }
    }
    assert!(receiver_gain > 0.0, "no cell received the film's mass");
    let before = species.cell(receiver)[ct::IN2];
    let after = advanced_species.cell(receiver)[ct::IN2];
    assert!(
        after > before,
        "the nitrogen the film brought did not land ({before} became {after})"
    );
    let total: f64 = advanced_species.cell(receiver).iter().sum();
    assert!(
        (total - 1.0).abs() < 1.0e-12,
        "the receiving cell's composition sums to {total}, not one"
    );

    // open: the case schema stores the coolant's choice as a string
    // and the cli turns it into this mixture with
    // film_cooling::coolant_composition, so the mapping asserted above
    // is the public part of that contract.
}

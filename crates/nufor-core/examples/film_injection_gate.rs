//! gate the blown-in film's injection accounting.
//!
//! the film COOLANT CROSSES the wall into the gas, so it carries mass,
//! energy and composition into the chamber, and it does so whatever the
//! wall itself is doing thermally. the gate measures three things:
//!
//! - the chamber's mass gain under a wall that states no temperature,
//!   which must be POSITIVE: a wall without a thermal treatment must
//!   not swallow the coolant;
//!
//! - the chamber's mass gain against the film's own book-keeping,
//!   which must MATCH: the report's coolant mass flow is the mass the
//!   chamber physically took;
//!
//! - the composition the film brought, which must LAND as nitrogen in
//!   the receiving cell.
//!
//! the numbers are printed as one json object so the python gate can
//! read them without parsing prose.

use nufor_core::combustion_march::{h2o2_ic_energy, CombustionState};
use nufor_core::combustion_thermo as ct;
use nufor_core::combustion_viscous_axi::{
    add_axi_combustion_viscous, axi_cell_volume, axi_combustion_transport, wall_diagnostics,
    AxiWallOpts,
};
use nufor_core::ConservedState2d;
use nufor_core::Grid2d;

/// the chamber the gate measures: eight by eight, with a flat wall at
/// y = 0.02 and the gas above it.
fn chamber() -> Grid2d {
    nufor_core::grid2d(8, 8, 0.0, 0.08, 0.0, 0.08).expect("the chamber's grid")
}

/// the resting gas: stoichiometric hydrogen and oxygen at seven
/// megapascals.
fn resting_gas(g: &Grid2d) -> (ConservedState2d, CombustionState) {
    let n = g.nx * g.ny;
    let rho = 12.64f64;
    let p = 7.0e6f64;
    let e_int = h2o2_ic_energy(rho, p);
    let state = ConservedState2d {
        rho: vec![rho; n],
        mx: vec![0.0; n],
        my: vec![0.0; n],
        e: vec![rho * e_int; n],
    };
    (state, CombustionState::stoich_h2o2(n))
}

/// the coolant's composition: pure nitrogen, at the closure's nitrogen
/// slot.
fn nitrogen() -> [f64; ct::NS] {
    let mut composition = [0.0f64; ct::NS];
    let slot = ct::SPECIES_NAMES
        .iter()
        .position(|&name| name == "N2")
        .expect("the closure has a nitrogen");
    composition[slot] = 1.0;
    composition
}

fn main() {
    let g = chamber();
    let nx = g.nx;
    let ny = g.ny;
    let n = nx * ny;
    let (state, species) = resting_gas(&g);
    let dist = |_x: f64, y: f64| y - 0.02;
    let dt = 1.0e-3;
    let mass_flux = 5.0f64;

    // the ADIABATIC wall: it states no temperature, so it does no
    // thermal exchange at all. the film still crosses it.
    let opts = AxiWallOpts {
        temperature: 0.0,
        emissivity: 0.0,
        film_coefficient: 0.0,
        film_injection_mass_flux: mass_flux,
        film_injection_temperature: 350.0,
        film_injection_composition: nitrogen(),
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

    // the chamber's mass gain, summed over the axisymmetric volumes.
    let gained: f64 = (0..n)
        .map(|k| (advanced.rho[k] - state.rho[k]) * axi_cell_volume(&g, k % nx, k / nx))
        .sum();

    // the film's own book-keeping.
    let report =
        wall_diagnostics(&state, &species, &g, 4, Some(&dist), opts).expect("the wall report");
    let booked = report.film_injection_mass_flow_kg_per_second * dt;

    // the nitrogen the film brought, as a mass fraction of the cell
    // that received the most.
    let nitrogen_slot = ct::SPECIES_NAMES
        .iter()
        .position(|&name| name == "N2")
        .expect("the closure has a nitrogen");
    let mut best_cell = 0usize;
    let mut best_gain = 0.0f64;
    for k in 0..n {
        let gain = (advanced.rho[k] - state.rho[k]) * axi_cell_volume(&g, k % nx, k / nx);
        if gain > best_gain {
            best_gain = gain;
            best_cell = k;
        }
    }
    let received_nitrogen = advanced_species.y[nitrogen_slot][best_cell];

    println!(
        "{{\"mass_flux_kg_per_m2_s\":{mass_flux:?},\"adiabatic_wall_temperature_k\":{:?},\
\"chamber_mass_gain_kg\":{gained:?},\"reported_mass_flow_kg_per_s\":{:?},\
\"booked_mass_kg\":{booked:?},\"relative_mass_error\":{:?},\
\"film_area_m2\":{:?},\"received_nitrogen_fraction\":{received_nitrogen:?}}}",
        opts.temperature,
        report.film_injection_mass_flow_kg_per_second,
        ((gained - booked) / booked).abs(),
        report.film_injection_area_square_metres,
    );
}

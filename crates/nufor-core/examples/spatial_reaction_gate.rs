//! gate finite-rate GRI-Mech source terms on a spatial composition.
//!
//! the gate lays a small field of parcels down, advances each cell with
//! the SPATIAL SOURCE PASS, and advances the same parcels with the
//! ZERO-DIMENSIONAL REACTOR MARCH. the two must agree cell by cell,
//! because the spatial pass is the reactor's own implicit step applied
//! per cell — anything else would mean the field's chemistry diverges
//! from the parcel reactor everyone calibrates against.
//!
//! the numbers print as one json object for the python gate.

use nufor_core::chem_react::apply_gri30_source;
use nufor_core::chem_spatial::SpatialComposition;
use nufor_core::gri30;
use nufor_core::gri30_reactor::{
    backward_euler_step, internal_energy_mass, temperature_at_energy, NS,
};

/// one parcel: a lean methane-air mixture at the given temperature.
fn parcel(dilution: f64, temperature: f64) -> ([f64; NS], f64) {
    let air = nufor_core::gri30_reactor::stoich_air("CH4", 1.0).expect("stoichiometric air");
    let ch4 = gri30::species_index("CH4").unwrap();
    let mut y = [0.0f64; NS];
    for s in 0..NS {
        y[s] = if s == ch4 {
            dilution
        } else {
            (1.0 - dilution) * air[s]
        };
    }
    let total: f64 = y.iter().sum();
    for value in y.iter_mut() {
        *value /= total;
    }
    let energy = internal_energy_mass(&y, temperature);
    (y, energy)
}

fn main() {
    // four parcels spanning cold-to-hot and rich-to-lean.
    let seeds = [
        (0.055, 900.0f64),
        (0.055, 1400.0),
        (0.09, 1200.0),
        (0.03, 1700.0),
    ];
    let ncells = seeds.len();
    let mut initial = Vec::with_capacity(ncells);
    let mut energies = Vec::with_capacity(ncells);
    for &(dilution, temperature) in &seeds {
        let (y, energy) = parcel(dilution, temperature);
        initial.push(y);
        energies.push(energy);
    }
    // the parcels themselves are indexed directly; no field matrix is
    // needed because the comparison is per cell.
    let densities = vec![1.0f64; ncells];

    // for each cell, find a step size the reactor's implicit solve
    // ACCEPTS, walking down from the nominal dt, then apply THAT SAME
    // dt ONCE on both the spatial pass and the reactor's own step. one
    // accepted step, identical on both sides, is the honest comparison.
    let nominal_dt = 1.0e-6f64;
    let mut total_reacted = 0usize;
    let mut total_deferred = 0usize;
    let mut worst_iterations = 0u32;
    let mut worst_relative = 0.0f64;
    let mut worst_abs = 0.0f64;
    let mut max_sum_error = 0.0f64;
    let mut chosen_dts = Vec::with_capacity(ncells);
    for cell in 0..ncells {
        let mut dt = nominal_dt;
        let mut accepted = None;
        for _ in 0..40 {
            let mut iterations = 0u32;
            match backward_euler_step(
                &initial[cell],
                densities[cell],
                energies[cell],
                dt,
                &mut iterations,
            ) {
                Ok(Some(stepped)) => {
                    accepted = Some((stepped, dt, iterations));
                    break;
                }
                Ok(None) => dt *= 0.5,
                Err(message) => panic!("the reactor's step errored on parcel {cell}: {message}"),
            }
        }
        let (expected, dt, iterations) =
            accepted.expect("the reactor's step never converged on any parcel");
        worst_iterations = worst_iterations.max(iterations);
        chosen_dts.push(dt);

        // the spatial pass on a one-cell field holding this parcel.
        let mut one_cell = SpatialComposition {
            y: (0..NS).map(|s| vec![initial[cell][s]]).collect(),
        };
        let report = apply_gri30_source(&mut one_cell, &[densities[cell]], &[energies[cell]], dt)
            .expect("source pass");
        total_reacted += report.cells_reacted;
        total_deferred += report.cells_deferred;
        assert_eq!(report.cells_reacted, 1, "cell {cell} deferred at dt {dt:e}");
        max_sum_error = max_sum_error.max(one_cell.max_sum_error());

        for (s, &want_raw) in expected.iter().enumerate() {
            let got = one_cell.y[s][0];
            let want = want_raw.max(0.0);
            let scale = want.abs().max(1.0e-12);
            let rel = (got - want).abs() / scale;
            if rel > worst_relative {
                worst_relative = rel;
            }
            worst_abs = worst_abs.max((got - want).abs());
        }
    }

    let temperatures_before: Vec<f64> = initial
        .iter()
        .zip(energies.iter())
        .map(|(parcel_, energy)| temperature_at_energy(parcel_, *energy))
        .collect();

    println!(
        "{{\"cells\":{ncells},\"cells_reacted\":{},\"cells_deferred\":{},\
\"worst_iterations\":{},\"max_sum_error\":{:?},\"worst_relative_error_vs_reactor\":{:?},\
\"chosen_dts_s\":{:?},\"worst_absolute_error_vs_reactor\":{:?},\"temperatures_before_k\":{:?}}}",
        total_reacted,
        total_deferred,
        worst_iterations,
        max_sum_error,
        worst_relative,
        chosen_dts,
        worst_abs,
        temperatures_before,
    );
}

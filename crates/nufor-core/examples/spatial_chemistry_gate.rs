//! gated measurement of the scheme-width spatial composition seam.
//!
//! the example drives the same public composition and transport path
//! used by the core tests. it measures a real GRI-Mech 3.0
//! composition crossing a planar flow step and writes the measurements
//! as json on standard output.

use nufor_core::chem_scheme::{ChemScheme, Gri30};
use nufor_core::chem_spatial::{advance_composition, SpatialComposition};

fn main() {
    let nx = 12usize;
    let ny = 6usize;
    let ncells = nx * ny;
    let methane = nufor_core::gri30::species_index("CH4").expect("CH4 is in GRI-Mech 3.0");
    let oxygen = nufor_core::gri30::species_index("O2").expect("O2 is in GRI-Mech 3.0");
    let nitrogen = nufor_core::gri30::species_index("N2").expect("N2 is in GRI-Mech 3.0");

    let mut fractions = vec![0.0; Gri30.n_species()];
    fractions[nitrogen] = 0.5;
    let background_species = Gri30.n_species() - 3;
    for (species, fraction) in fractions.iter_mut().enumerate() {
        if species != methane && species != oxygen && species != nitrogen {
            *fraction = 0.5 / background_species as f64;
        }
    }
    let mut composition = SpatialComposition::uniform(&Gri30, ncells, &fractions)
        .expect("the GRI-Mech composition sums to one");

    for (species, row) in composition.y.iter_mut().enumerate() {
        for (cell, value) in row.iter_mut().enumerate() {
            let inlet = cell % nx == 0;
            let methane_fraction = if inlet { 0.20 } else { 0.005 };
            let oxygen_fraction = if inlet { 0.35 } else { 0.195 };
            let nitrogen_fraction = 1.0 - methane_fraction - oxygen_fraction;
            *value = if species == methane {
                methane_fraction
            } else if species == oxygen {
                oxygen_fraction
            } else if species == nitrogen {
                nitrogen_fraction * 0.5
            } else {
                nitrogen_fraction * 0.5 / background_species as f64
            };
        }
    }

    let inlet_methane = composition.y[methane][0];
    let receiver_cell = 1;
    let initial_receiver_methane = composition.y[methane][receiver_cell];

    let density = vec![1.0; ncells];
    let dxs = vec![0.01; nx];
    let dys = vec![0.01; ny];
    let mut fx = vec![vec![0.0; nx + 1]; ny];
    let fy = vec![vec![0.0; ny + 1]; nx];
    for row in &mut fx {
        for face in row.iter_mut() {
            *face = 0.2;
        }
    }

    let residual =
        advance_composition(&mut composition, &density, 0.05, &dxs, &dys, &fx, &fy, None)
            .expect("the spatial composition step succeeds");
    let received_methane = composition.y[methane][receiver_cell];
    let minimum_positive_fraction = composition
        .y
        .iter()
        .flatten()
        .fold(f64::INFINITY, |minimum, value| minimum.min(*value));

    println!(
        "{{\"scheme\":\"{}\",\"species\":{},\"cells\":{},\"inlet_methane_fraction\":{:.12e},\
\"initial_receiver_methane_fraction\":{:.12e},\"received_methane_fraction\":{:.12e},\
\"transferred_methane_fraction\":{:.12e},\"maximum_composition_sum_error\":{:.12e},\
\"minimum_positive_species_fraction\":{:.12e},\"transport_relative_update\":{:.12e}}}",
        Gri30.name(),
        Gri30.n_species(),
        ncells,
        inlet_methane,
        initial_receiver_methane,
        received_methane,
        received_methane - initial_receiver_methane,
        composition.max_sum_error(),
        minimum_positive_fraction,
        residual,
    );
}

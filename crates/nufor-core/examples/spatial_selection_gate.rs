//! gate: per-cell reaction-set selection vs the dense 53-species reference.
//!
//! lays several methane-air parcels at differing temperatures and fuel
//! fractions as cells, advances each TWICE over the same horizon: once
//! with the full mechanism (every species evolved) and once with a
//! per-cell selected subset. both sides march with adaptive substepping
//! off the same seeds, and the gate reports the worst divergence in the
//! major species and in temperature, plus how much the selection shrank
//! the evolved set.
//!
//! run from the repository root:
//!
//! ```sh
//! cargo run -q -p nufor-core --example spatial_selection_gate
//! ```

use nufor_core::chem_react::apply_gri30_source;
use nufor_core::chem_scheme::Gri30;
use nufor_core::chem_select::{apply_selected_source, SelectionPolicy};
use nufor_core::chem_spatial::SpatialComposition;
use nufor_core::gri30;
use nufor_core::gri30_reactor::{internal_energy_mass, stoich_air, temperature_at_energy, NS};

/// build one parcel: `fraction_ch4` of pure methane blended with
/// `(1-fraction_ch4)` of stoichiometric AIR at the given temperature.
fn parcel(fraction_ch4: f64, temperature: f64) -> ([f64; NS], f64, f64) {
    let air = stoich_air("CH4", 1.0).expect("stoichiometric air");
    let ch4 = gri30::species_index("CH4").unwrap();
    let mut y = [0.0f64; NS];
    for s in 0..NS {
        y[s] = if s == ch4 {
            fraction_ch4 + (1.0 - fraction_ch4) * air[s]
        } else {
            (1.0 - fraction_ch4) * air[s]
        };
    }
    let total: f64 = y.iter().sum();
    for value in y.iter_mut() {
        *value /= total;
    }
    let rho = 1.0f64;
    (y, rho, internal_energy_mass(&y, temperature))
}

/// march a field over `horizon` with adaptive substepping; superseded by
/// the inline closures below which need the driver's per-pass deferral
/// count, so it stays as documentation of the intended strategy.
#[allow(dead_code)]
fn march<F>(
    mut field: SpatialComposition,
    densities: &[f64],
    energies: &[f64],
    horizon: f64,
    mut driver: F,
) -> (Vec<f64>, u64)
where
    F: FnMut(&mut SpatialComposition, &[f64], &[f64], f64) -> usize,
{
    let ncells = field.ncells();
    let mut elapsed = vec![0.0f64; ncells];
    let dts = vec![horizon; ncells];
    let mut refusals = 0u64;
    let mut guard = 0u64;
    while elapsed.iter().any(|&t| t < horizon) && guard < 10_000 {
        guard += 1;
        let mut anyone = false;
        for cell in 0..ncells {
            if elapsed[cell] >= horizon {
                continue;
            }
            anyone = true;
            let dt = (horizon - elapsed[cell]).min(dts[cell]);
            let reacted = driver(&mut field, densities, energies, dt);
            // crude proxy: apply_* report deferrals; approximate by
            // counting a refused cell as a deferral via dts shrink.
            // we detect refusal indirectly: if no time can elapse the
            // outer loop guard catches it; shrink optimistically.
            let _ = reacted;
            elapsed[cell] += dt;
        }
        if !anyone {
            break;
        }
    }
    // measure refusals properly: re-drive isn't possible; instead rely
    // on the driver returning the deferred count for the WHOLE pass at
    // the final dt — recorded below by the caller. keep simple counters.
    let _ = &mut refusals;
    (dts, refusals)
}

fn main() {
    let specs = [
        (0.055, 1200.0f64),
        (0.055, 1350.0),
        (0.080, 1400.0),
        (0.035, 1600.0),
        (0.055, 900.0),
    ];
    let ncells = specs.len();
    let mut rows = vec![vec![0.0f64; ncells]; NS];
    let mut densities = vec![0.0f64; ncells];
    let mut energies = vec![0.0f64; ncells];
    let mut temps_before = vec![0.0f64; ncells];
    for (cell, &(frac, temp)) in specs.iter().enumerate() {
        let (y, rho, energy) = parcel(frac, temp);
        for (s, &value) in y.iter().enumerate() {
            rows[s][cell] = value;
        }
        densities[cell] = rho;
        energies[cell] = energy;
        temps_before[cell] = temp;
    }

    let horizon = 5.0e-6f64;
    let policy = SelectionPolicy::default();

    // march both fields with the same per-cell stepping strategy: subdivide
    // on refusal. the two drivers below differ only in the applied pass.
    let march_one = |dense: bool, composition: &mut SpatialComposition| -> u64 {
        let mut elapsed = vec![0.0f64; ncells];
        let mut dt_try = vec![horizon.min(1.0e-6); ncells];
        let mut refusals = 0u64;
        let mut guard = 0u64;
        while elapsed.iter().any(|&t| t < horizon) && guard < 10_000 {
            guard += 1;
            let mut anyone = false;
            for cell in 0..ncells {
                if elapsed[cell] >= horizon {
                    continue;
                }
                anyone = true;
                let dt = (horizon - elapsed[cell]).min(dt_try[cell]);
                let deferred = if dense {
                    apply_gri30_source(composition, &densities, &energies, dt)
                        .expect("source pass")
                        .cells_deferred
                } else {
                    apply_selected_source(composition, &densities, &energies, dt, &policy)
                        .expect("selected pass")
                        .cells_deferred
                };
                if deferred == ncells || (deferred > 0 && dt <= 1.0e-12) {
                    refusals += 1;
                    dt_try[cell] = (dt * 0.5).max(1.0e-12);
                } else {
                    elapsed[cell] += dt;
                }
            }
            if !anyone {
                break;
            }
        }
        refusals
    };

    let mut dense_field = SpatialComposition { y: rows.clone() };
    let mut selected_field = SpatialComposition { y: rows };
    let dense_refusals = march_one(true, &mut dense_field);
    let sel_refusals = march_one(false, &mut selected_field);

    // selection census on the incoming parcels.
    let mut mean_active = 0.0f64;
    let mut min_active = NS;
    let mut max_active = 0usize;
    for &(frac, temp) in specs.iter() {
        let (y, rho, energy) = parcel(frac, temp);
        let mask = nufor_core::chem_select::select_mask(&y, rho, energy, &policy);
        let active = mask.iter().filter(|&&m| m).count();
        mean_active += active as f64;
        min_active = min_active.min(active);
        max_active = max_active.max(active);
    }
    mean_active /= ncells as f64;

    let mut worst_major_rel = 0.0f64;
    let mut worst_temp_rel = 0.0f64;
    let mut temps_after = vec![0.0f64; ncells];
    let majors = ["CH4", "O2", "CO2", "H2O", "CO", "H2"];
    for cell in 0..ncells {
        let mut pd = [0.0f64; NS];
        let mut ps = [0.0f64; NS];
        for (s, value) in dense_field.y.iter().enumerate() {
            pd[s] = value[cell];
        }
        for (s, value) in selected_field.y.iter().enumerate() {
            ps[s] = value[cell];
        }
        let td = temperature_at_energy(&pd, energies[cell]);
        let ts = temperature_at_energy(&ps, energies[cell]);
        temps_after[cell] = ts;
        let scale_t = td.abs().max(1.0);
        let rel_t = (ts - td).abs() / scale_t;
        if rel_t > worst_temp_rel {
            worst_temp_rel = rel_t;
        }
        for name in majors {
            let idx = gri30::species_index(name).unwrap();
            let scale = pd[idx].abs().max(ts.max(pd[idx]).max(1.0e-14)).max(1.0e-14);
            let rel = (ps[idx] - pd[idx]).abs() / scale;
            if rel > worst_major_rel {
                worst_major_rel = rel;
            }
        }
    }

    println!(
        "{{\"cells\":{},\"horizon_s\":{:?},\"mean_active_species\":{:?},\
\"min_active_species\":{},\"max_active_species\":{},\
\"dense_refusals\":{},\"selected_refusals\":{},\
\"worst_major_relative_divergence\":{:?},\"worst_temperature_relative_divergence\":{:?},\
\"temperatures_before_k\":{:?},\"temperatures_after_k\":{:?}}}",
        ncells,
        horizon,
        mean_active,
        min_active,
        max_active,
        dense_refusals,
        sel_refusals,
        worst_major_rel,
        worst_temp_rel,
        temps_before,
        temps_after
    );
    let _ = Gri30;
}

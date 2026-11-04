//! gate: the threaded slab sharder partitions cells and matches the serial
//! accumulation, exercised on a wide chemical source pass.
//!
//! demonstrates the multi-core primitive the chamber march will adopt:
//! per-cell chemistry applied over disjoint slabs reproduces the serial
//! result bit-for-bit because each cell's update depends only on its own
//! state.
//!
//! run from the repository root:
//!
//! ```sh
//! cargo run -q -p nufor-core --example par_slabs_gate
//! ```

use nufor_core::chem_scheme::Gri30;
use nufor_core::chem_select::{apply_selected_source, SelectionPolicy};
use nufor_core::chem_spatial::SpatialComposition;
use nufor_core::gri30;
use nufor_core::gri30_reactor::{internal_energy_mass, stoich_air, NS};
use nufor_core::par_slabs::run_slabs;

fn parcel(fraction_ch4: f64, temperature: f64) -> ([f64; NS], f64, f64) {
    let air = stoich_air("CH4", 1.0).unwrap();
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
    for v in y.iter_mut() {
        *v /= total;
    }
    (y, 1.0, internal_energy_mass(&y, temperature))
}

/// one slab's completed chemistry pass: the cell range it covered, the
/// advanced composition rows, and how many of its cells were deferred.
struct SlabResult {
    start: usize,
    end: usize,
    y: Vec<Vec<f64>>,
    deferred: usize,
}

fn main() {
    // a field of 257 cells at varied states: odd size exercises ragged slabs.
    let ncells = 257usize;
    let mut rows = vec![vec![0.0f64; ncells]; NS];
    let mut densities = vec![1.0f64; ncells];
    let mut energies = vec![0.0f64; ncells];
    for cell in 0..ncells {
        let frac = 0.03 + 0.05 * ((cell as f64 * 0.37) % 1.0);
        let temp = 1100.0 + 500.0 * ((cell as f64 * 0.61) % 1.0);
        let (y, rho, energy) = parcel(frac, temp);
        for (s, &v) in y.iter().enumerate() {
            rows[s][cell] = v;
        }
        densities[cell] = rho;
        energies[cell] = energy;
    }

    let dt = 5.0e-7f64;
    let policy = SelectionPolicy::default();

    // serial reference.
    let mut serial = SpatialComposition { y: rows.clone() };
    let serial_report =
        apply_selected_source(&mut serial, &densities, &energies, dt, &policy).unwrap();

    // threaded: partition the field into slabs and apply the selected pass
    // on each DISJOINT sub-field independently, reassembling afterwards.
    // valid because per-cell chemistry depends only on the cell's own
    // composition, so the union of slab results equals the serial result
    // bit for bit. dispatched through the pool so the slabs really run in
    // parallel.
    let bounds: Vec<(usize, usize)> = {
        let threads = 4;
        let chunk = ncells.div_ceil(threads);
        let mut b = Vec::new();
        let mut start = 0;
        while start < ncells {
            let end = (start + chunk).min(ncells);
            b.push((start, end));
            start = end;
        }
        b
    };

    // run each slab's pass inside the pool, capturing its sub-field result.
    // run_slabs hands each worker a contiguous range of SLAB indices, so the
    // closure maps that range to one output per slab.
    let slab_outputs: Vec<SlabResult> = {
        let bounds_ref = &bounds;
        let rows_ref = &rows;
        let dens_ref = &densities;
        let ener_ref = &energies;
        run_slabs(bounds.len(), 4, move |a, b| {
            (a..b)
                .flat_map(|bi| {
                    let (start, end) = bounds_ref[bi];
                    let mut sub = SpatialComposition {
                        y: (0..NS).map(|s| rows_ref[s][start..end].to_vec()).collect(),
                    };
                    let rep = apply_selected_source(
                        &mut sub,
                        &dens_ref[start..end],
                        &ener_ref[start..end],
                        dt,
                        &policy,
                    )
                    .unwrap();
                    vec![SlabResult {
                        start,
                        end,
                        y: sub.y,
                        deferred: rep.cells_deferred,
                    }]
                })
                .collect::<Vec<_>>()
        })
        .into_iter()
        .flatten()
        .collect()
    };

    let mut assembled = vec![vec![0.0f64; ncells]; NS];
    let mut deferred_total = 0usize;
    for slab in slab_outputs {
        let SlabResult {
            start,
            end,
            y,
            deferred,
        } = slab;
        deferred_total += deferred;
        for (s, row) in y.iter().enumerate() {
            assembled[s][start..end].copy_from_slice(row);
        }
    }
    let threaded = SpatialComposition { y: assembled };

    let mut worst_bit_diff = 0.0f64;
    let mut differs = 0usize;
    for s in 0..NS {
        for k in 0..ncells {
            let a = serial.y[s][k];
            let b = threaded.y[s][k];
            if a.to_bits() != b.to_bits() {
                differs += 1;
                worst_bit_diff = worst_bit_diff.max((a - b).abs());
            }
        }
    }
    let sum_drift = threaded.max_sum_error();

    println!(
        "{{\"cells\":{},\"slabs\":{},\"serial_deferred\":{},\"threaded_deferred\":{},\
\"bitwise_equal_entries\":{},\"differing_entries\":{},\"worst_abs_diff_when_differing\":{:?},\
\"threaded_max_sum_error\":{:?}}}",
        ncells,
        bounds.len(),
        serial_report.cells_deferred,
        deferred_total,
        ncells * NS - differs,
        differs,
        worst_bit_diff,
        sum_drift
    );
    let _ = Gri30;
}

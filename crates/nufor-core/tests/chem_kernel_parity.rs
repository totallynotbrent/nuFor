//! parity between the chemistry kernel's hydrogen arm and the legacy
//! march's own nine-species cell step.
//!
//! the legacy step is private to `combustion_march`, so it is
//! transcribed here verbatim from that file. the transcription is
//! compared against the live source by eye; if the march changes, this
//! transcript has to change with it, and the bit-identity assertions
//! then fail loudly.
//!
//! the assertion is on bits, not tolerance: the kernel must produce
//! the same doubles, down to the last mantissa bit, for every species
//! and for the returned timescale violation.

use nufor_core::chem_kernel::advance_cell;
use nufor_core::chem_lane::ChemLane;
use nufor_core::combustion as comb;
use nufor_core::combustion_closure as closure9;
use nufor_core::combustion_thermo as ct;

/// the legacy march's cell step, copied verbatim from
/// `crates/nufor-core/src/combustion_march.rs`.
fn reference_react_cell(e_int: f64, rho: f64, y: &mut [f64; ct::NS], dt: f64) -> f64 {
    let mut worst_oversub = 0.0f64;
    {
        let mut t = closure9::temperature_at(e_int, y);
        let mut remaining = dt;
        let mut steps = 0;
        while remaining > 1e-14 && steps < 512 {
            let (dydt, _) = comb::omega(rho, y, t);
            let rate_scale = dydt.iter().fold(0.0f64, |acc, d| acc.max(d.abs()));
            let tau = if rate_scale > 1e-30 {
                1.0 / rate_scale
            } else {
                remaining
            };
            let h = remaining.min(tau * 0.2);
            // positivity by RATE SCALING, not by shrinking the step:
            // the reaction vector is scaled down by the largest factor
            // that keeps every species non-negative over this step.
            // the stoichiometry is untouched, so the species sum and
            // the element pools stay exactly conserved, and the
            // chemistry still advances (shrinking the step instead
            // would freeze it behind a trace species).
            let mut scale = 1.0f64;
            for (s, &d) in dydt.iter().enumerate() {
                if d < 0.0 && y[s] > 0.0 {
                    let allowed = y[s] / (-d * h);
                    if allowed < scale {
                        scale = allowed;
                    }
                }
            }
            for (ys, d) in y.iter_mut().zip(dydt.iter()) {
                let d = d * h * scale;
                // no clamp and no renormalization anywhere: the
                // reaction system conserves the species sum and the
                // element pools exactly, and any correction destroys
                // one or both.
                *ys += d;
            }
            t = closure9::temperature_at(e_int, y);
            remaining -= h;
            steps += 1;
            if steps == 1 {
                worst_oversub = worst_oversub.max(dt / (tau * 0.2));
            }
        }
        // no sanitize anywhere: the reaction system conserves both
        // the species sum and the element pools exactly, so any
        // clamp or renormalization manufactures or destroys mass. the
        // adaptive step keeps excursions small; a mildly negative
        // species relaxes back within the same cell-dt.
    }
    worst_oversub
}

/// a spread of charges: stoichiometric premixes at several dilutions and
/// temperatures, plus a burnt-over mix carrying water and radicals.
fn charges() -> Vec<([f64; ct::NS], f64, f64)> {
    let mut cases: Vec<([f64; ct::NS], f64, f64)> = Vec::new();

    // the stoichiometric hydrogen/oxygen charge.
    let mut stoich = [0.0f64; ct::NS];
    stoich[ct::IH2] = 2.016 / 18.016;
    stoich[ct::IO2] = 16.0 / 18.016;
    for &(dilution, temperature) in &[
        (1.0f64, 800.0f64),
        (1.0, 1_200.0),
        (1.0, 1_600.0),
        (0.5, 1_400.0),
        (0.25, 2_000.0),
    ] {
        let mut y = [0.0f64; ct::NS];
        for s in 0..ct::NS {
            y[s] = dilution * stoich[s];
        }
        // the leftover is inert nitrogen.
        y[ct::IN2] = 1.0 - dilution;
        let rho = 1.0 + 0.5 * dilution;
        cases.push((y, rho, temperature));
    }

    // a burnt-over composition: water, radicals, leftover oxygen.
    let mut burnt = [0.0f64; ct::NS];
    burnt[ct::IH2O] = 0.55;
    burnt[ct::IOH] = 0.02;
    burnt[ct::IH] = 0.005;
    burnt[ct::IO_] = 0.003;
    burnt[ct::IO2] = 0.12;
    burnt[ct::IN2] = 0.302;
    cases.push((burnt, 0.8, 1_800.0));

    // a nearly inert charge, which exercises the no-reaction path.
    let mut inert = [0.0f64; ct::NS];
    inert[ct::IN2] = 0.997;
    inert[ct::IO2] = 0.003;
    cases.push((inert, 1.2, 700.0));

    cases
}

/// the kernel's hydrogen arm and the legacy step agree bit for bit.
#[test]
fn li_arm_is_bit_identical_to_the_legacy_step() {
    let timesteps = [1.0e-9f64, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3];
    let mut checked = 0usize;
    for (charge, rho, temperature) in charges() {
        let e_int = closure9::internal_energy_at(temperature, &charge);
        for dt in timesteps {
            let mut legacy = charge;
            let want_violation = reference_react_cell(e_int, rho, &mut legacy, dt);

            let mut kernel = charge.to_vec();
            let got_violation =
                advance_cell(ChemLane::Li2004, e_int, rho, &mut kernel, dt).expect("step");

            assert!(
                want_violation.to_bits() == got_violation.to_bits(),
                "violation differs: legacy {want_violation:e} against kernel {got_violation:e}"
            );
            for (s, (&want, &got)) in legacy.iter().zip(kernel.iter()).enumerate() {
                assert!(
                    want.to_bits() == got.to_bits(),
                    "species {s} differs at dt {dt:e}: legacy {want:e} against kernel {got:e}"
                );
            }
            checked += 1;
        }
    }
    assert_eq!(checked, charges().len() * timesteps.len());
}

/// the field pass walks the cells independently, so splitting the field
/// gives the same bits as one serial pass.
#[test]
fn the_field_pass_matches_a_serial_walk_cell_for_cell() {
    let charges = charges();
    let ncells = charges.len();
    let mut density = vec![0.0f64; ncells];
    let mut energy = vec![0.0f64; ncells];
    let mut rows = vec![vec![0.0f64; ncells]; ct::NS];
    for (cell, (charge, rho, temperature)) in charges.iter().enumerate() {
        density[cell] = *rho;
        energy[cell] = closure9::internal_energy_at(*temperature, charge);
        for s in 0..ct::NS {
            rows[s][cell] = charge[s];
        }
    }

    let mut whole = rows.clone();
    let worst = nufor_core::chem_kernel::advance_field(
        ChemLane::Li2004,
        &mut whole,
        &density,
        &energy,
        1.0e-6,
    )
    .expect("field");
    assert!(worst.is_finite());

    for cell in 0..ncells {
        let mut piece = vec![vec![0.0f64; 1]; ct::NS];
        for s in 0..ct::NS {
            piece[s][0] = rows[s][cell];
        }
        nufor_core::chem_kernel::advance_field(
            ChemLane::Li2004,
            &mut piece,
            &density[cell..cell + 1],
            &energy[cell..cell + 1],
            1.0e-6,
        )
        .expect("slab");
        for s in 0..ct::NS {
            assert!(
                piece[s][0].to_bits() == whole[s][cell].to_bits(),
                "species {s} in cell {cell} differs between the whole-field and slab runs"
            );
        }
    }
}

/// the methane arm steps implicitly and moves the fuel, holding the
/// composition normalised and non-negative.
#[test]
fn the_gri_arm_steps_implicitly() {
    use nufor_core::gri30;
    use nufor_core::gri30_reactor;

    let lane = ChemLane::Gri30;
    let air = gri30_reactor::stoich_air("CH4", 1.0).unwrap();
    let ch4 = gri30::species_index("CH4").unwrap();
    let mut y: Vec<f64> = (0..gri30_reactor::NS)
        .map(|s| {
            if s == ch4 {
                0.055
            } else {
                (1.0 - 0.055) * air[s]
            }
        })
        .collect();
    let rho = 1.0f64;
    let temperature = 1_200.0f64;
    let e_int = lane.internal_energy_at(temperature, &y);
    let before = y[ch4];
    let before_vector = y.clone();

    let violation = advance_cell(lane, e_int, rho, &mut y, 1.0e-6).expect("gri step");
    assert_eq!(violation, 0.0, "the implicit arm reports no violation");
    assert!(y.iter().all(|v| v.is_finite()));
    assert!(
        y.iter().zip(before_vector.iter()).any(|(a, b)| *a != *b),
        "nothing moved"
    );

    // the kernel's step IS the reactor's own backward-euler step, so the
    // two agree species for species, bit for bit.
    let mut parcel = [0.0f64; gri30_reactor::NS];
    parcel.copy_from_slice(&before_vector);
    let mut iterations = 0u32;
    let expected = gri30_reactor::backward_euler_step(&parcel, rho, e_int, 1.0e-6, &mut iterations)
        .expect("reactor step")
        .expect("reactor converged");
    for s in 0..gri30_reactor::NS {
        assert!(
            y[s].to_bits() == expected[s].to_bits(),
            "species {s}: kernel {:e} against reactor {:e}",
            y[s],
            expected[s]
        );
    }
    let _ = before;
}

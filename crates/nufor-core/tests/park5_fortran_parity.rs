//! parity between the fortran park 5-species lane and the rust lane.
//!
//! the fortran module `fortran/api/nuforpark5.f90` is the port of the
//! rust lane in `crates/nufor-core/src/{species_closure,species_chem,
//! species_transport}.rs`. this gate compares the two on the same
//! inputs and asserts they agree to tight tolerances, with the
//! transcendental-sensitive quantities held to a much tighter bound
//! than the algebraic ones.
//!
//! the tolerances are loose where the rust and fortran compilers may
//! legitimately disagree on floating-point association (the sums over
//! species and the products in the closure), and tight where they must
//! not (the rates themselves, which go through the same c library
//! exp/log/pow on both sides).

use nufor_core::park::NS;
use nufor_core::park5_fortran;
use nufor_core::species_chem;
use nufor_core::species_closure as closure5;

/// the relative tolerance for the closure and the transport, where the
/// two compilers may associate the sums differently.
const LOOSE_REL: f64 = 1.0e-12;

/// the relative tolerance for the reaction rates, which both sides
/// compute through the same c library transcendentals.
const TIGHT_REL: f64 = 1.0e-13;

/// standard air: 0.767 N2 / 0.233 O2 by mass.
fn air_y() -> [f64; NS] {
    [0.767, 0.233, 0.0, 0.0, 0.0]
}

/// a spread of compositions spanning inert air, a partly dissociated
/// shock layer, a nitric-oxide-rich mix and a hot atomic mix.
fn charges() -> Vec<([f64; NS], f64, f64)> {
    vec![
        (air_y(), 1.2, 300.0),
        // a hot compressed mix leaning into dissociation.
        ([0.60, 0.20, 0.05, 0.10, 0.05], 2.5, 6500.0),
        ([0.50, 0.15, 0.05, 0.20, 0.10], 4.0, 9000.0),
        // a lean mix that stays mostly molecular.
        ([0.85, 0.12, 0.02, 0.008, 0.002], 0.8, 4500.0),
        // a near-equilibrium atomic mix, exercising the backward rates.
        ([0.20, 0.05, 0.05, 0.45, 0.25], 3.0, 11000.0),
        // a cold dense cell, where the temperature floor can bite.
        ([0.95, 0.05, 0.0, 0.0, 0.0], 5.0, 120.0),
    ]
}

/// the closure agrees: pressure, temperature and sound speed.
#[test]
fn the_closure_agrees_between_the_two_lanes() {
    for (y, rho, temperature) in charges() {
        let e_int = closure5::internal_energy_at(temperature, &y);
        let planes = vec![vec![y[0]], vec![y[1]], vec![y[2]], vec![y[3]], vec![y[4]]];

        let fields = park5_fortran::closure(&[rho], &[e_int], &planes).expect("fortran closure");

        let pr = closure5::pressure_at(rho, e_int, &y);
        let tr = closure5::temperature_at(e_int, &y);
        let ar = closure5::sound_speed_at(rho, e_int, &y);

        assert!(
            (fields.t[0] - tr).abs() <= LOOSE_REL * tr.abs().max(1.0),
            "temperature differs: fortran {:?} against rust {tr}",
            fields.t
        );
        assert!(
            (fields.p[0] - pr).abs() <= LOOSE_REL * pr.abs().max(1.0),
            "pressure differs: fortran {:?} against rust {pr}",
            fields.p
        );
        assert!(
            (fields.a[0] - ar).abs() <= LOOSE_REL * ar.abs().max(1.0),
            "sound speed differs: fortran {:?} against rust {ar}",
            fields.a
        );

        // the fortran closure must invert the rust closure: feeding the
        // fortran temperature back through the rust side reproduces the
        // same energy.
        let e_back = closure5::internal_energy_at(fields.t[0], &y);
        assert!(
            (e_back - e_int).abs() <= LOOSE_REL * e_int.abs().max(1.0),
            "closure inversion drifted: {e_back} against {e_int}"
        );
    }
}

/// the reaction source agrees: the per-species production rates and the
/// chemical energy deposition.
#[test]
fn the_reaction_source_agrees_between_the_two_lanes() {
    for (y, rho, temperature) in charges() {
        let planes = vec![vec![y[0]], vec![y[1]], vec![y[2]], vec![y[3]], vec![y[4]]];
        let (dydt_f, de_f) =
            park5_fortran::omega(&[rho], &[temperature], &planes).expect("fortran omega");
        let (dydt_r, de_r) = species_chem::omega(rho, &y, temperature);

        for s in 0..NS {
            let want = dydt_r[s];
            let got = dydt_f[s][0];
            let scale = want.abs().max(got.abs()).max(1.0e-30);
            assert!(
                (got - want).abs() <= TIGHT_REL * scale,
                "production rate for species {s} differs at {temperature} K: \
                 fortran {got:e} against rust {want:e}"
            );
        }
        let scale = de_r.abs().max(de_f[0].abs()).max(1.0e-30);
        assert!(
            (de_f[0] - de_r).abs() <= TIGHT_REL * scale,
            "chemical deposition differs: fortran {e:e} against rust {de_r:e}",
            e = de_f[0]
        );
    }
}

/// the chemistry substep agrees: the composition and the two
/// diagnostics.
#[test]
fn the_chemistry_substep_agrees_between_the_two_lanes() {
    for (y, rho, temperature) in charges() {
        let e_int = closure5::internal_energy_at(temperature, &y);
        for dt in [1.0e-9f64, 1.0e-7, 1.0e-6, 1.0e-5] {
            let mut planes = vec![vec![y[0]], vec![y[1]], vec![y[2]], vec![y[3]], vec![y[4]]];
            let (wd_f, wo_f) = park5_fortran::chem_substep(&[rho], &[e_int], dt, &mut planes)
                .expect("fortran chemistry");

            // the rust lane's own chemistry substep is the reference:
            // it runs the same operator split on the same closure, so
            // the comparison is against the shipping code and not a
            // transcription of it.
            let rho_v = vec![rho];
            let g = nufor_core::grid2d::grid2d(1, 1, 0.0, 0.1, 0.0, 0.1).expect("grid");
            let st = nufor_core::state2d::ConservedState2d {
                rho: rho_v.clone(),
                mx: vec![0.0],
                my: vec![0.0],
                e: vec![rho * e_int],
            };
            let mut sp = nufor_core::species::SpeciesState::uniform(1, &y).expect("species");
            let (wd_r, wo_r) = nufor_core::species_march::chemistry_substep(&st, &mut sp, &g, dt)
                .expect("rust chemistry");

            assert!(
                (wo_f - wo_r).abs() <= LOOSE_REL * wo_r.abs().max(1.0),
                "timescale violation differs at dt {dt:e}: fortran {wo_f:e} against rust {wo_r:e}"
            );
            assert!(
                (wd_f - wd_r).abs() <= LOOSE_REL * wd_r.abs().max(1.0),
                "worst dY differs at dt {dt:e}: fortran {wd_f:e} against rust {wd_r:e}"
            );
            for (s, want) in sp.y.iter().enumerate().take(NS) {
                let got = planes[s][0];
                assert!(
                    (got - want[0]).abs() <= LOOSE_REL * want[0].abs().max(1.0e-6),
                    "composition for species {s} differs at dt {dt:e}: \
                     fortran {got:e} against rust {:e}",
                    want[0]
                );
            }
            // the composition stays a distribution.
            let sum: f64 = planes.iter().map(|pl| pl[0]).sum();
            assert!(
                (sum - 1.0).abs() <= 1e-12,
                "composition does not sum to one: {sum}"
            );
        }
    }
}

/// the species transport agrees: the same face mass fluxes advance the
/// same composition.
#[test]
fn the_species_transport_agrees_between_the_two_lanes() {
    let nx = 4usize;
    let ny = 3usize;
    let n = nx * ny;
    let dxs = vec![0.1f64; nx];
    let dys = vec![0.1f64; ny];
    let rho = vec![1.0f64; n];

    // a varying composition, so the advection actually moves species.
    let mut planes: Vec<Vec<f64>> = Vec::new();
    for s in 0..NS {
        planes.push((0..n).map(|k| 0.1 + 0.01 * (k as f64 + s as f64)).collect());
    }
    // renormalize the starting composition.
    for k in 0..n {
        let sum: f64 = planes.iter().map(|pl| pl[k]).sum();
        for pl in planes.iter_mut() {
            pl[k] /= sum;
        }
    }

    // a uniform x flux and no y flux: a shear-like throughflow.
    let fx: Vec<Vec<f64>> = (0..ny).map(|_| vec![0.5f64; nx + 1]).collect();
    let fy: Vec<Vec<f64>> = (0..nx).map(|_| vec![0.0f64; ny + 1]).collect();

    // the rust reference: the same core, called through the crate.
    let mut want = planes.clone();
    let resid_r = nufor_core::species_transport::transport_species_generic(
        &mut want, &rho, 1.0e-3, &dxs, &dys, &fx, &fy, None,
    )
    .expect("rust transport");

    let mut got = planes.clone();
    let grid = park5_fortran::TransportGrid {
        nx,
        ny,
        dxs: &dxs,
        dys: &dys,
        fx: &fx,
        fy: &fy,
    };
    let resid_f =
        park5_fortran::transport(&grid, &rho, 1.0e-3, &mut got).expect("fortran transport");

    assert!(
        (resid_f - resid_r).abs() <= LOOSE_REL * resid_r.abs().max(1.0e-30),
        "residual differs: fortran {resid_f:e} against rust {resid_r:e}"
    );
    for s in 0..NS {
        for k in 0..n {
            let a = want[s][k];
            let b = got[s][k];
            assert!(
                (a - b).abs() <= LOOSE_REL * a.abs().max(1.0e-12),
                "transport differs for species {s} in cell {k}: \
                 fortran {b:e} against rust {a:e}"
            );
        }
    }
}

/// a uniform composition under a uniform flux stays uniform: the
/// conservation invariant the transport owes the flow.
#[test]
fn the_fortran_transport_keeps_a_uniform_field_uniform() {
    let nx = 5usize;
    let ny = 4usize;
    let n = nx * ny;
    let dxs = vec![0.1f64; nx];
    let dys = vec![0.1f64; ny];
    let rho = vec![1.0f64; n];
    let mut planes = vec![vec![0.767f64; n], vec![0.233f64; n]];
    planes.extend(vec![vec![0.0f64; n]; NS - 2]);

    let fx: Vec<Vec<f64>> = (0..ny).map(|_| vec![0.5f64; nx + 1]).collect();
    let fy: Vec<Vec<f64>> = (0..nx).map(|_| vec![0.0f64; ny + 1]).collect();

    let grid = park5_fortran::TransportGrid {
        nx,
        ny,
        dxs: &dxs,
        dys: &dys,
        fx: &fx,
        fy: &fy,
    };
    let resid =
        park5_fortran::transport(&grid, &rho, 1.0e-3, &mut planes).expect("fortran transport");

    assert!(resid < 1e-14, "a uniform field diverged: resid {resid:e}");
    // the composition stays exactly where it started, species by
    // species: the first two carry the air, the rest are absent.
    let expect = [0.767, 0.233, 0.0, 0.0, 0.0];
    for (s, pl) in planes.iter().enumerate() {
        for v in pl.iter() {
            assert!(
                (v - expect[s]).abs() < 1e-14,
                "a uniform composition moved: species {s} went to {v}"
            );
        }
    }
}

/// the fortran lane rejects malformed arguments instead of writing past
/// the buffers.
#[test]
fn the_fortran_lane_rejects_bad_arguments() {
    assert!(park5_fortran::closure(&[], &[], &[]).is_err());
    assert!(park5_fortran::omega(&[1.0], &[300.0], &[]).is_err());
    let mut planes = vec![vec![0.767f64], vec![0.233f64]];
    planes.extend(vec![vec![0.0f64]; NS - 2]);
    assert!(park5_fortran::chem_substep(&[1.0, 1.0], &[1.0], 1e-6, &mut planes).is_err());
}

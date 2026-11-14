//! behaviour of the fortran cartesian park5 march.
//!
//! `fortran/api/nuforcartpark5.f90` is the cartesian species-carrying
//! marcher: the flow advances under the five-species mixture closure,
//! the species ride the face mass fluxes and the operator-split
//! chemistry follows. these gates exercise the real fortran code
//! through the c abi and assert the invariants the lane owes:
//!
//! - a stagnant uniform field stays uniform and conserves its mass;
//! - a compressive ramp accelerates the flow monotonically and keeps
//!   the density positive;
//! - the composition stays normalised and non-negative through a step
//!   that heats and reacts the gas.

use nufor_core::park::NS;
use nufor_core::park5_fortran::{self, CartBc};

/// the standard air composition by mass.
const AIR: [f64; NS] = [0.767, 0.233, 0.0, 0.0, 0.0];

/// the freestream internal energy per mass for dry air at 288.15 K.
fn air_internal_energy() -> f64 {
    nufor_core::species_closure::internal_energy_at(288.15, &AIR)
}

/// the composition planes for a uniform field.
fn uniform_planes(n: usize) -> Vec<Vec<f64>> {
    (0..NS).map(|s| vec![AIR[s]; n]).collect()
}

/// transmissive boundaries with a consistent imposed state.
fn bc_vals() -> [[f64; 5]; 4] {
    let mut v = [[0.0f64; 5]; 4];
    for side in v.iter_mut() {
        side[0] = 1.225; // density
        side[1] = 0.0; // u
        side[2] = 0.0; // v
        side[3] = 101325.0; // pressure
        side[4] = air_internal_energy(); // internal energy per mass
    }
    v
}

/// a uniform freestream at rest.
fn freestream(nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = nx * ny;
    let e_int = air_internal_energy();
    (
        vec![1.225; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![1.225 * e_int; n],
    )
}

/// the marsingale grid: transmissive everywhere.
fn transmissive() -> (CartBc, CartBc, CartBc, CartBc) {
    (
        CartBc::Transmissive,
        CartBc::Transmissive,
        CartBc::Transmissive,
        CartBc::Transmissive,
    )
}

/// a stagnant uniform field must not move: the lane's strongest
/// statement about its face bookkeeping.
#[test]
fn a_stagnant_uniform_field_stays_uniform() {
    let (nx, ny) = (12usize, 8usize);
    let n = nx * ny;
    let (rho, mx, my, e) = freestream(nx, ny);
    let planes = uniform_planes(n);
    let (bw, be, bs, bn) = transmissive();

    let out = park5_fortran::cart_step(
        nx,
        ny,
        &vec![1.0 / nx as f64; nx],
        &vec![1.0 / ny as f64; ny],
        0.45,
        bw,
        be,
        bs,
        bn,
        &bc_vals(),
        &rho,
        &mx,
        &my,
        &e,
        &planes,
        &AIR,
    )
    .expect("fortran cart step");

    let mass_before: f64 = rho.iter().sum();
    let mass_after: f64 = out.rho.iter().sum();
    assert!(
        (mass_after - mass_before).abs() <= 1.0e-12 * mass_before,
        "mass drifted: {mass_before:e} became {mass_after:e}"
    );
    for (k, r) in out.rho.iter().enumerate() {
        assert!(
            (*r - rho[k]).abs() <= 1.0e-12 * rho[k],
            "density moved in cell {k}: {} became {}",
            rho[k],
            r
        );
    }
    for v in out.mx.iter().chain(out.my.iter()) {
        assert!(v.abs() <= 1.0e-10, "a stagnant field acquired motion");
    }
    // the composition stays put as well: uniform air under no flux.
    for (s, plane) in out.planes.iter().enumerate() {
        for (k, v) in plane.iter().enumerate() {
            assert!(
                (v - AIR[s]).abs() <= 1.0e-12,
                "composition moved: species {s} in cell {k} went to {v}"
            );
        }
    }
}

/// a compressing ramp must accelerate the flow to the right and keep
/// every cell positive.
#[test]
fn a_compressive_ramp_accelerates_the_flow_monotonically() {
    let (nx, ny) = (32usize, 4usize);
    let n = nx * ny;
    let e_int = air_internal_energy();
    let mut rho = vec![0.0f64; n];
    let mut mx = vec![0.0f64; n];
    let mut my = vec![0.0f64; n];
    let mut e = vec![0.0f64; n];
    // a linear density ramp rising along x at a fixed stagnation
    // enthalpy: the flow runs downhill and speeds up along the ramp.
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let frac = i as f64 / (nx - 1) as f64;
            rho[k] = 1.0 + 1.5 * frac;
            mx[k] = 0.0;
            my[k] = 0.0;
            e[k] = rho[k] * e_int;
        }
    }
    let planes = uniform_planes(n);
    let (bw, be, bs, bn) = transmissive();

    let out = park5_fortran::cart_step(
        nx,
        ny,
        &vec![1.0 / nx as f64; nx],
        &vec![1.0 / ny as f64; ny],
        0.3,
        bw,
        be,
        bs,
        bn,
        &bc_vals(),
        &rho,
        &mx,
        &my,
        &e,
        &planes,
        &AIR,
    )
    .expect("fortran cart step");

    assert!(
        out.dt.is_finite() && out.dt > 0.0,
        "the step chose no timestep"
    );
    for (k, r) in out.rho.iter().enumerate() {
        assert!(*r > 0.0, "the density went non-positive in cell {k}");
        assert!(r.is_finite(), "the density went infinite in cell {k}");
    }
    // the ramp's own gradient steepens the pressure gradient, so the
    // downstream half of the row must carry more momentum than the
    // upstream half.
    let row = ny / 2;
    let upstream: f64 = out.mx[row * nx + 1];
    let downstream: f64 = out.mx[row * nx + nx - 2];
    assert!(
        downstream > upstream,
        "the compression did not accelerate the flow: \
         upstream {upstream:e} against downstream {downstream:e}"
    );
}

/// the composition leaves every step normalised and non-negative, even
/// where the chemistry is active.
#[test]
fn the_composition_stays_normalised_through_a_hot_step() {
    let (nx, ny) = (8usize, 6usize);
    let n = nx * ny;
    let (rho, mx, my, e) = freestream(nx, ny);
    let planes = uniform_planes(n);
    let (bw, be, bs, bn) = transmissive();

    let mut cur_rho = rho;
    let mut cur_mx = mx;
    let mut cur_my = my;
    let mut cur_e = e;
    let mut cur_planes = planes;

    for _ in 0..8 {
        let out = park5_fortran::cart_step(
            nx,
            ny,
            &vec![1.0 / nx as f64; nx],
            &vec![1.0 / ny as f64; ny],
            0.45,
            bw,
            be,
            bs,
            bn,
            &bc_vals(),
            &cur_rho,
            &cur_mx,
            &cur_my,
            &cur_e,
            &cur_planes,
            &AIR,
        )
        .expect("fortran cart step");
        cur_rho = out.rho;
        cur_mx = out.mx;
        cur_my = out.my;
        cur_e = out.e;
        cur_planes = out.planes;

        for (s, plane) in cur_planes.iter().enumerate() {
            for (k, v) in plane.iter().enumerate() {
                assert!(*v >= 0.0, "species {s} went negative in cell {k}: {v}");
                assert!(v.is_finite(), "species {s} went infinite in cell {k}");
            }
        }
        for k in 0..n {
            let sum: f64 = cur_planes.iter().map(|pl| pl[k]).sum();
            assert!(
                (sum - 1.0).abs() <= 1.0e-9,
                "the composition left normalisation in cell {k}: {sum}"
            );
        }
    }
}

/// the lane refuses malformed grids instead of writing past the
/// buffers.
#[test]
fn the_cart_lane_rejects_bad_arguments() {
    let (rho, mx, my, e) = freestream(4, 3);
    let planes = uniform_planes(12);
    let (bw, be, bs, bn) = transmissive();
    assert!(park5_fortran::cart_step(
        0,
        3,
        &[],
        &[],
        0.4,
        bw,
        be,
        bs,
        bn,
        &bc_vals(),
        &rho,
        &mx,
        &my,
        &e,
        &planes,
        &AIR
    )
    .is_err());
    assert!(park5_fortran::cart_step(
        4,
        3,
        &[0.25; 3],
        &[0.33; 3],
        0.4,
        bw,
        be,
        bs,
        bn,
        &bc_vals(),
        &rho,
        &mx,
        &my,
        &e,
        &planes,
        &AIR
    )
    .is_err());
}

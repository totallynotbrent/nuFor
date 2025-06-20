//! integration tests for CFL time-step control, pinned to the published formula.

use nufor_core::{cfl_dt, Error};

// reference max characteristic speed reimplemented from the definition.
fn reference_speed(gamma: f64, rho: f64, m: f64, e: f64) -> f64 {
    let u = m / rho;
    let p = (gamma - 1.0) * (e - 0.5 * rho * u * u);
    let a = (gamma * p / rho).sqrt();
    u.abs() + a
}

#[test]
fn static_gas_gives_the_sound_speed_bound() {
    // gamma = 2, rho = 1, p = 8: a = 4, u = 0, so s_max = 4 and
    // dt = 0.5 * 1 / 4 = 0.125, all exact in binary.
    let step = cfl_dt(2.0, 0.5, 1.0, &[1.0], &[0.0], &[8.0]).unwrap();
    assert!((step.max_speed - 4.0).abs() < 1e-12);
    assert!((step.dt - 0.125).abs() < 1e-12);
}

#[test]
fn moving_gas_adds_the_flow_speed() {
    // rho = 4, u = 2, p = 8 with gamma = 2: a = 2, u = 2, so s_max = 4.
    let step = cfl_dt(2.0, 0.5, 1.0, &[4.0], &[8.0], &[16.0]).unwrap();
    assert!((step.max_speed - 4.0).abs() < 1e-12);
    assert!((step.dt - 0.125).abs() < 1e-12);
}

#[test]
fn supersonic_state_matches_the_reference_formula() {
    // u = 5, rho = 1, p = 1 with gamma = 1.4: s_max = 5 + sqrt(1.4).
    let step = cfl_dt(1.4, 0.9, 2.0, &[1.0], &[5.0], &[15.0]).unwrap();
    let expected_speed = reference_speed(1.4, 1.0, 5.0, 15.0);
    assert!((step.max_speed - expected_speed).abs() < 1e-12);
    assert!((step.dt - 0.9 * 2.0 / expected_speed).abs() < 1e-12);
}

#[test]
fn the_fastest_cell_sets_the_global_step() {
    // two cells: a quiet one and a fast one; the global step follows the fast cell.
    let fast_speed = reference_speed(1.4, 1.0, 5.0, 15.0);
    let step = cfl_dt(1.4, 0.5, 1.0, &[1.0, 1.0], &[0.0, 5.0], &[2.5, 15.0]).unwrap();
    assert!((step.max_speed - fast_speed).abs() < 1e-12);
    assert!((step.dt - 0.5 / fast_speed).abs() < 1e-12);
}

#[test]
fn dt_scales_with_cfl_and_dx() {
    // halving the Courant number halves dt; doubling the cell width doubles it.
    let base = cfl_dt(2.0, 0.5, 1.0, &[1.0], &[0.0], &[8.0]).unwrap();
    let smaller = cfl_dt(2.0, 0.25, 1.0, &[1.0], &[0.0], &[8.0]).unwrap();
    let wider = cfl_dt(2.0, 0.5, 2.0, &[1.0], &[0.0], &[8.0]).unwrap();
    assert!((smaller.dt - 0.5 * base.dt).abs() < 1e-15);
    assert!((wider.dt - 2.0 * base.dt).abs() < 1e-15);
    assert!((smaller.max_speed - base.max_speed).abs() < 1e-15);
}

#[test]
fn cfl_dt_rejects_out_of_range_parameters() {
    let (rho, m, e) = ([1.0], [0.0], [2.5]);
    for cfl in [0.0, -0.1, 1.5, f64::NAN] {
        assert_eq!(
            cfl_dt(1.4, cfl, 1.0, &rho, &m, &e),
            Err(Error::InvalidArgs),
            "cfl = {cfl} should fail"
        );
    }
    for dx in [0.0, -1.0, f64::NAN] {
        assert_eq!(
            cfl_dt(1.4, 0.5, dx, &rho, &m, &e),
            Err(Error::InvalidArgs),
            "dx = {dx} should fail"
        );
    }
    for gamma in [1.0, 0.9, f64::NAN] {
        assert_eq!(
            cfl_dt(gamma, 0.5, 1.0, &rho, &m, &e),
            Err(Error::InvalidArgs),
            "gamma = {gamma} should fail"
        );
    }
}

#[test]
fn cfl_dt_rejects_nonphysical_states() {
    // non-positive density on either cell.
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[0.0], &[0.0], &[2.5]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[-1.0], &[0.0], &[2.5]),
        Err(Error::KernelFailure)
    );
    // u = 2 needs e_t > 2; at or below it has no positive internal energy.
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[1.0], &[2.0], &[2.0]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[1.0], &[2.0], &[1.9]),
        Err(Error::KernelFailure)
    );
}

#[test]
fn cfl_dt_rejects_non_finite_state_components() {
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[f64::NAN], &[0.0], &[2.5]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[1.0], &[f64::INFINITY], &[2.5]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[1.0], &[0.0], &[f64::NAN]),
        Err(Error::KernelFailure)
    );
}

#[test]
fn cfl_dt_rejects_mismatched_or_empty_slices() {
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[1.0], &[0.0, 1.0], &[2.5]),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &[1.0], &[0.0], &[2.5, 3.0]),
        Err(Error::InvalidArgs)
    );
    let empty: [f64; 0] = [];
    assert_eq!(
        cfl_dt(1.4, 0.5, 1.0, &empty, &empty, &empty),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn cfl_dt_handles_a_single_cell() {
    // A lone cell must give its own speed and step, not a constant.
    let expected_speed = reference_speed(1.4, 1.0, 5.0, 15.0);
    let step = cfl_dt(1.4, 0.7, 0.25, &[1.0], &[5.0], &[15.0]).unwrap();
    assert!((step.max_speed - expected_speed).abs() < 1e-12);
    assert!((step.dt - 0.7 * 0.25 / expected_speed).abs() < 1e-12);
}

//! 2d grid, state, and equation-of-state checks that extend the 1d ones.

use nufor_core::{
    check_physical2d, cons_to_prim2d, eos_mach2d, eos_pressure2d, eos_sound_speed2d, grid2d,
    prim_to_cons2d, Error, Grid2d,
};

#[test]
fn grid2d_geometry_is_uniform_and_closed() {
    let g: Grid2d = grid2d(4, 3, 0.0, 1.0, 0.0, 2.0).unwrap();
    assert_eq!(g.nx, 4);
    assert_eq!(g.ny, 3);
    assert_eq!(g.centers_x.len(), 12);
    assert!((g.dx - 0.25).abs() < 1e-13);
    assert!((g.dy - 2.0 / 3.0).abs() < 1e-13);
    // first cell is (0.125, 1/3); last is (0.875, 5/3).
    assert!((g.centers_x[0] - 0.125).abs() < 1e-13);
    assert!((g.centers_y[0] - 1.0 / 3.0).abs() < 1e-13);
    assert!((g.centers_x[11] - 0.875).abs() < 1e-13);
    assert!((g.centers_y[11] - 5.0 / 3.0).abs() < 1e-13);
    // faces close the domain exactly.
    assert_eq!(g.faces_x.len(), 5);
    assert_eq!(g.faces_y.len(), 4);
    assert!((g.faces_x[0] - 0.0).abs() < 1e-13 && (g.faces_x[4] - 1.0).abs() < 1e-13);
    assert!((g.faces_y[3] - 2.0).abs() < 1e-13);
}

#[test]
fn grid2d_rejects_degenerate_inputs() {
    assert_eq!(
        grid2d(0, 2, 0.0, 1.0, 0.0, 1.0).unwrap_err(),
        Error::InvalidArgs
    );
    assert_eq!(
        grid2d(2, 0, 0.0, 1.0, 0.0, 1.0).unwrap_err(),
        Error::InvalidArgs
    );
    assert_eq!(
        grid2d(2, 2, 1.0, 1.0, 0.0, 1.0).unwrap_err(),
        Error::InvalidArgs
    );
    assert_eq!(
        grid2d(2, 2, 0.0, 1.0, 2.0, 1.0).unwrap_err(),
        Error::InvalidArgs
    );
    assert_eq!(
        grid2d(2, 2, f64::NAN, 1.0, 0.0, 1.0).unwrap_err(),
        Error::InvalidArgs
    );
}

#[test]
fn prim_to_cons2d_uses_definition_values() {
    let (mx, my, e) = prim_to_cons2d(&[2.0, 3.0], &[3.0, 1.0], &[4.0, 2.0], &[8.0, 6.0]).unwrap();
    assert_eq!(mx, vec![6.0, 3.0]);
    assert_eq!(my, vec![8.0, 6.0]);
    assert_eq!(e, vec![16.0, 18.0]);
}

#[test]
fn cons_to_prim2d_round_trips() {
    let (u, v, et) = cons_to_prim2d(&[2.0, 3.0], &[6.0, 9.0], &[4.0, 6.0], &[16.0, 24.0]).unwrap();
    assert!((u[0] - 3.0).abs() < 1e-12 && (u[1] - 3.0).abs() < 1e-12);
    assert!((v[0] - 2.0).abs() < 1e-12 && (v[1] - 2.0).abs() < 1e-12);
    assert!((et[0] - 8.0).abs() < 1e-12 && (et[1] - 8.0).abs() < 1e-12);
}

#[test]
fn eos_pressure2d_uses_definition_values() {
    // p = 0.4 * 2 * (10 - 0.5*(2^2+2^2)) = 0.8 * 6 = 4.8
    let p = eos_pressure2d(1.4, &[2.0], &[10.0], &[2.0], &[2.0]).unwrap();
    assert!((p[0] - 4.8).abs() < 1e-12);
}

#[test]
fn eos_sound_speed_and_mach() {
    let a = eos_sound_speed2d(1.4, &[4.8], &[2.0]).unwrap();
    assert!((a[0] - (1.4f64 * 4.8 / 2.0).sqrt()).abs() < 1e-12);
    let m = eos_mach2d(&[3.0, 1.0], &[4.0, 0.0], &[5.0, 5.0]).unwrap();
    assert!((m[0] - 1.0).abs() < 1e-12);
    assert!((m[1] - 0.2).abs() < 1e-12);
}

#[test]
fn eos_pressure2d_rejects_bad_density() {
    assert_eq!(
        eos_pressure2d(1.4, &[-1.0], &[8.0], &[0.0], &[0.0]).unwrap_err(),
        Error::InvalidArgs
    );
    assert_eq!(
        eos_pressure2d(1.0, &[1.0], &[8.0], &[0.0], &[0.0]).unwrap_err(),
        Error::InvalidArgs
    );
}

#[test]
fn check_flags_bad_2d_states_and_passes_healthy_ones() {
    use nufor_core::ConservedState2d;
    let healthy = ConservedState2d {
        rho: vec![1.0, 2.0],
        mx: vec![0.5, 1.0],
        my: vec![0.0, 0.0],
        e: vec![2.5, 5.0],
    };
    assert!(check_physical2d(&healthy).ok);
    let mut bad = healthy.clone();
    bad.rho[1] = f64::NAN;
    let chk = check_physical2d(&bad);
    assert!(!chk.ok);
    assert_eq!(chk.bad_cell, Some(1));
    let mut neg = healthy.clone();
    neg.rho[0] = -0.25;
    assert!(!check_physical2d(&neg).ok);
}

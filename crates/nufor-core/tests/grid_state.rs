//! integration tests for uniform 1D grid geometry and the conserved/primitive conversion.

use nufor_core::{cons_to_prim, grid1d, prim_to_cons, Error};

#[test]
fn grid_is_uniform_and_closes_the_domain() {
    let g = grid1d(8, 0.0, 4.0).unwrap();
    assert_eq!(g.dx, 0.5);
    assert_eq!(g.centers.len(), 8);
    assert_eq!(g.faces.len(), 9);
    // centers sit at (i - 1/2) * dx
    let expected: Vec<f64> = (0..8).map(|i| (i as f64 + 0.5) * 0.5).collect();
    assert_eq!(g.centers, expected);
    // faces close the domain exactly, interior faces at multiples of dx
    assert_eq!(g.faces.first(), Some(&0.0));
    assert_eq!(g.faces.last(), Some(&4.0));
    assert_eq!(g.faces[3], 1.5);
    assert_eq!(g.faces[7], 3.5);
}

#[test]
fn grid_offsets_shift_with_xmin() {
    let g = grid1d(4, 2.0, 3.0).unwrap();
    assert_eq!(g.dx, 0.25);
    assert_eq!(g.centers[0], 2.125);
    assert_eq!(g.faces.first(), Some(&2.0));
    assert_eq!(g.faces.last(), Some(&3.0));
}

#[test]
fn grid_rejects_degenerate_input() {
    assert_eq!(grid1d(1, 0.0, 1.0), Err(Error::InvalidArgs));
    assert_eq!(grid1d(0, 0.0, 1.0), Err(Error::InvalidArgs));
    assert_eq!(grid1d(4, 1.0, 1.0), Err(Error::InvalidArgs));
    assert_eq!(grid1d(4, 2.0, 1.0), Err(Error::InvalidArgs));
    assert!(grid1d(4, f64::NAN, 1.0).is_err());
}

#[test]
fn prim_to_cons_matches_the_definitions() {
    // powers of two keep every value exact in binary.
    let rho = [4.0, 2.0, 8.0];
    let u = [2.0, -3.0, 0.5];
    let et = [5.0, 7.0, 10.0];
    let (m, e) = prim_to_cons(&rho, &u, &et).unwrap();
    assert_eq!(m, vec![8.0, -6.0, 4.0]);
    assert_eq!(e, vec![20.0, 14.0, 80.0]);
}

#[test]
fn cons_to_prim_inverts_prim_to_cons() {
    let rho = [4.0, 2.0, 8.0];
    let u = [2.0, -3.0, 0.5];
    let et = [5.0, 7.0, 10.0];
    let (m, e) = prim_to_cons(&rho, &u, &et).unwrap();
    let (u2, et2) = cons_to_prim(&rho, &m, &e).unwrap();
    assert_eq!(u2, u.to_vec());
    assert_eq!(et2, et.to_vec());
}

#[test]
fn round_trip_recovers_primitives_within_machine_epsilon() {
    // the identity only holds for binary-representable values; u = m/rho then rho*u is not bit-exact.
    let rho = [1.4, 0.9, 3.3];
    let u = [0.7, -1.1, 2.9];
    let et = [11.3, 4.7, 15.1];
    let (m, e) = prim_to_cons(&rho, &u, &et).unwrap();
    let (u2, et2) = cons_to_prim(&rho, &m, &e).unwrap();
    for (a, b) in u2.iter().zip(u.iter()) {
        assert!((a - b).abs() < 1e-12);
    }
    for (a, b) in et2.iter().zip(et.iter()) {
        assert!((a - b).abs() < 1e-12);
    }
}

#[test]
fn negative_or_zero_density_is_a_numerical_failure() {
    let bad_rho = [4.0, -1.0, 8.0];
    let zero_rho = [0.0, 2.0, 8.0];
    let u = [2.0, 2.0, 2.0];
    let et = [5.0, 5.0, 5.0];
    assert_eq!(prim_to_cons(&bad_rho, &u, &et), Err(Error::KernelFailure));
    assert_eq!(prim_to_cons(&zero_rho, &u, &et), Err(Error::KernelFailure));
    let m = [8.0, 8.0, 8.0];
    let e = [20.0, 20.0, 20.0];
    assert_eq!(cons_to_prim(&bad_rho, &m, &e), Err(Error::KernelFailure));
}

#[test]
fn conversion_rejects_mismatched_or_empty_slices() {
    let rho = [1.0, 2.0];
    let u = [1.0];
    let et = [3.0, 4.0];
    assert_eq!(prim_to_cons(&rho, &u, &et), Err(Error::InvalidArgs));

    let rho_single = [1.0];
    let m = [2.0, 3.0];
    let e = [5.0];
    assert_eq!(cons_to_prim(&rho_single, &m, &e), Err(Error::InvalidArgs));

    let empty: [f64; 0] = [];
    assert_eq!(
        prim_to_cons(&empty, &empty, &empty),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn conversions_handle_single_cell() {
    let (m, e) = prim_to_cons(&[3.0], &[2.0], &[9.0]).unwrap();
    assert_eq!(m, vec![6.0]);
    assert_eq!(e, vec![27.0]);
    let (u, et) = cons_to_prim(&[3.0], &m, &e).unwrap();
    assert_eq!(u, vec![2.0]);
    assert_eq!(et, vec![9.0]);
}

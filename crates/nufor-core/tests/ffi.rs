//! FFI integration tests (spec 39): the Rust wrapper against the real Fortran
//! kernel library, linked by build.rs through CMake.

use nufor_core::{saxpy, version, Error};

#[test]
fn version_is_semver_string() {
    // Parity with the version string baked into nuforkernels.f90.
    assert_eq!(version(), "0.1.0");
}

#[test]
fn saxpy_matches_reference_loop() {
    let x = [1.0, 2.0, 3.0, 4.0];
    let mut y = [10.0, 20.0, 30.0, 40.0];
    let expected: Vec<f64> = y.iter().zip(x.iter()).map(|(y, x)| 2.5 * x + y).collect();

    saxpy(2.5, &x, &mut y).unwrap();

    assert_eq!(y.to_vec(), expected);
}

#[test]
fn saxpy_rejects_mismatched_lengths() {
    let x = [1.0, 2.0];
    let mut y = [1.0, 2.0, 3.0];
    assert_eq!(saxpy(1.0, &x, &mut y), Err(Error::InvalidArgs));
}

#[test]
fn saxpy_on_empty_slices_is_a_no_op() {
    let x: [f64; 0] = [];
    let mut y: [f64; 0] = [];
    saxpy(1.0, &x, &mut y).unwrap(); // must not call the kernel
}

#[test]
fn saxpy_ok_on_single_element() {
    let x = [7.0];
    let mut y = [1.0];
    saxpy(1.0, &x, &mut y).unwrap();
    assert_eq!(y[0], 8.0);
}

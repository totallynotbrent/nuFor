//! Step-6 FFI tests: HLL numerical flux for the 1D Euler equations (spec 52,
//! 90) with Davis wave-speed estimates. Expected values are either exact
//! physical fluxes (supersonic / identical-state cases) or the published HLL
//! formula applied to the input states (subsonic cases).

use nufor_core::{hll_flux, Error, HllFlux};

// Reference HLL flux for a face, reimplemented from the published formula so
// the kernel is pinned to the source rather than to a magic constant.
fn reference_hll(
    gamma: f64,
    rl: f64,
    ml: f64,
    el: f64,
    rr: f64,
    mr: f64,
    er: f64,
) -> (f64, f64, f64) {
    let ul = ml / rl;
    let pl = (gamma - 1.0) * (el - 0.5 * rl * ul * ul);
    let al = (gamma * pl / rl).sqrt();
    let ur = mr / rr;
    let pr = (gamma - 1.0) * (er - 0.5 * rr * ur * ur);
    let ar = (gamma * pr / rr).sqrt();
    let sl = (ul - al).min(ur - ar);
    let sr = (ul + al).max(ur + ar);
    let fl = (rl * ul, rl * ul * ul + pl, ul * (el + pl));
    let fr = (rr * ur, rr * ur * ur + pr, ur * (er + pr));
    if sl >= 0.0 {
        fl
    } else if sr <= 0.0 {
        fr
    } else {
        let w = 1.0 / (sr - sl);
        (
            w * (sr * fl.0 - sl * fr.0 + sl * sr * (rr - rl)),
            w * (sr * fl.1 - sl * fr.1 + sl * sr * (mr - ml)),
            w * (sr * fl.2 - sl * fr.2 + sl * sr * (er - el)),
        )
    }
}

fn assert_flux(flux: &HllFlux, face: usize, expected: (f64, f64, f64)) {
    assert!((flux.rho[face] - expected.0).abs() < 1e-12);
    assert!((flux.m[face] - expected.1).abs() < 1e-12);
    assert!((flux.e[face] - expected.2).abs() < 1e-12);
}

#[test]
fn identical_states_return_the_physical_flux() {
    // gamma = 2 and p = 8 make every value exact in binary.
    let flux = hll_flux(2.0, &[1.0], &[0.0], &[8.0], &[1.0], &[0.0], &[8.0]).unwrap();
    assert_flux(&flux, 0, (0.0, 8.0, 0.0));
    // A moving identical state: rho=4, u=2, p=8 with gamma=2 has e_t = 4.
    let flux = hll_flux(2.0, &[4.0], &[8.0], &[16.0], &[4.0], &[8.0], &[16.0]).unwrap();
    assert_flux(&flux, 0, (8.0, 24.0, 48.0));
}

#[test]
fn supersonic_flow_uses_the_left_physical_flux() {
    // u = 5, rho = 1, p = 1 with gamma = 1.4: u - a > 0, so both wave speeds
    // are positive and the flux is the left physical flux (exact).
    let flux = hll_flux(1.4, &[1.0], &[5.0], &[15.0], &[0.5], &[2.0], &[8.0]).unwrap();
    assert_flux(&flux, 0, (5.0, 26.0, 80.0));
}

#[test]
fn supersonic_flow_uses_the_right_physical_flux() {
    // u = -5, rho = 1, p = 1 with gamma = 1.4: u + a < 0, so both wave speeds
    // are negative and the flux is the right physical flux (exact).
    let flux = hll_flux(1.4, &[0.5], &[-2.0], &[8.0], &[1.0], &[-5.0], &[15.0]).unwrap();
    assert_flux(&flux, 0, (-5.0, 26.0, -80.0));
}

#[test]
fn hll_recovers_the_sod_face_flux() {
    // Sod face at t = 0: left (rho=1, u=0, p=1), right (rho=0.125, u=0,
    // p=0.1), gamma = 1.4. The momentum flux lands on 0.55 exactly because
    // the Davis speeds are symmetric about zero.
    let flux = hll_flux(1.4, &[1.0], &[0.0], &[2.5], &[0.125], &[0.0], &[0.25]).unwrap();
    assert_flux(
        &flux,
        0,
        reference_hll(1.4, 1.0, 0.0, 2.5, 0.125, 0.0, 0.25),
    );
    assert!((flux.m[0] - 0.55).abs() < 1e-12);
}

#[test]
fn subsonic_faces_match_the_reference_formula() {
    // A few generic subsonic faces and a low-density state.
    let cases = [
        (1.4, 1.0, 1.0, 3.0, 0.5, -0.25, 3.28125),
        (1.4, 2.0, -2.0, 1.875, 1.2, 0.48, 2.35),
        (1.4, 0.125, 0.05, 0.265625, 1.0, -1.5, 3.375),
        (1.67, 0.8, 0.4, 2.0, 1.2, -0.6, 3.5),
    ];
    for (gamma, rl, ml, el, rr, mr, er) in cases {
        let flux = hll_flux(gamma, &[rl], &[ml], &[el], &[rr], &[mr], &[er]).unwrap();
        assert_flux(&flux, 0, reference_hll(gamma, rl, ml, el, rr, mr, er));
    }
}

#[test]
fn many_faces_are_processed_in_one_call() {
    // Face 1 is Sod (subsonic middle state), face 2 is left-supersonic, and
    // face 3 is right-supersonic; each output slot must hold its own flux.
    let flux = hll_flux(
        1.4,
        &[1.0, 1.0, 0.5],
        &[0.0, 5.0, -2.0],
        &[2.5, 15.0, 8.0],
        &[0.125, 0.5, 1.0],
        &[0.0, 2.0, -5.0],
        &[0.25, 8.0, 15.0],
    )
    .unwrap();
    assert_flux(
        &flux,
        0,
        reference_hll(1.4, 1.0, 0.0, 2.5, 0.125, 0.0, 0.25),
    );
    assert_flux(&flux, 1, (5.0, 26.0, 80.0));
    assert_flux(&flux, 2, (-5.0, 26.0, -80.0));
}

#[test]
fn hll_rejects_nonpositive_density_on_either_side() {
    let (rho_r, m_r, e_r) = ([0.125], [0.0], [0.25]);
    assert_eq!(
        hll_flux(1.4, &[0.0], &[0.0], &[2.5], &rho_r, &m_r, &e_r),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        hll_flux(1.4, &[-1.0], &[0.0], &[2.5], &rho_r, &m_r, &e_r),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        hll_flux(1.4, &[1.0], &[0.0], &[2.5], &[0.0], &[0.0], &[0.25]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        hll_flux(1.4, &[1.0], &[0.0], &[2.5], &[-2.0], &[0.0], &[0.25]),
        Err(Error::KernelFailure)
    );
}

#[test]
fn hll_rejects_nonphysical_internal_energy() {
    // u = 2 needs e_t > 2 for positive internal energy; at or below it is
    // not an admissible Euler state on either side of the face.
    assert_eq!(
        hll_flux(1.4, &[1.0], &[2.0], &[2.0], &[1.0], &[0.0], &[2.5]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        hll_flux(1.4, &[1.0], &[2.0], &[1.9], &[1.0], &[0.0], &[2.5]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        hll_flux(1.4, &[1.0], &[0.0], &[2.5], &[1.0], &[2.0], &[2.0]),
        Err(Error::KernelFailure)
    );
}

#[test]
fn hll_rejects_non_finite_state_components() {
    assert_eq!(
        hll_flux(1.4, &[f64::NAN], &[0.0], &[2.5], &[1.0], &[0.0], &[2.5]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        hll_flux(
            1.4,
            &[1.0],
            &[f64::INFINITY],
            &[2.5],
            &[1.0],
            &[0.0],
            &[2.5]
        ),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        hll_flux(1.4, &[1.0], &[0.0], &[2.5], &[1.0], &[0.0], &[f64::NAN]),
        Err(Error::KernelFailure)
    );
}

#[test]
fn hll_rejects_invalid_gamma() {
    let (rho, m, e) = ([1.0], [0.0], [2.5]);
    assert_eq!(
        hll_flux(1.0, &rho, &m, &e, &rho, &m, &e),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        hll_flux(0.9, &rho, &m, &e, &rho, &m, &e),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        hll_flux(f64::NAN, &rho, &m, &e, &rho, &m, &e),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn hll_rejects_mismatched_or_empty_slices() {
    assert_eq!(
        hll_flux(
            1.4,
            &[1.0],
            &[0.0],
            &[2.5],
            &[0.125, 0.125],
            &[0.0],
            &[0.25],
        ),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        hll_flux(1.4, &[1.0], &[0.0, 1.0], &[2.5], &[0.125], &[0.0], &[0.25]),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        hll_flux(1.4, &[1.0], &[0.0], &[2.5], &[0.125], &[0.0], &[0.25, 0.5]),
        Err(Error::InvalidArgs)
    );
    let empty: [f64; 0] = [];
    assert_eq!(
        hll_flux(1.4, &empty, &empty, &empty, &empty, &empty, &empty),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn hll_handles_a_single_face() {
    // A lone subsonic face must match the reference formula, not a constant.
    let flux = hll_flux(2.0, &[2.0], &[2.0], &[9.0], &[1.0], &[2.0], &[5.0]).unwrap();
    assert_flux(&flux, 0, reference_hll(2.0, 2.0, 2.0, 9.0, 1.0, 2.0, 5.0));
}

//! Step-5 FFI tests: ideal-gas equation of state (pressure, sound speed,
//! Mach, temperature) and its physical-validity checks (spec 25, 94, 138).

use nufor_core::{
    cons_to_prim, eos_mach, eos_pressure, eos_sound_speed, eos_temperature, prim_to_cons, Error,
};

#[test]
fn pressure_matches_the_definition() {
    // gamma = 2 makes gamma - 1 = 1 so every value below is exact in binary.
    let rho = [4.0, 2.0, 8.0];
    let u = [2.0, -3.0, 0.5];
    let et = [5.0, 7.0, 10.0];
    let p = eos_pressure(2.0, &rho, &et, &u).unwrap();
    // e_int = e_t - u^2/2 -> 3, 2.5, 9.875; p = (gamma-1)*rho*e_int
    assert_eq!(p, vec![12.0, 5.0, 79.0]);
}

#[test]
fn pressure_recovers_the_sod_left_state() {
    // Sod left state rho=1, u=0, p=1 with gamma=1.4 has e_t = p/(rho*(gamma-1)).
    let p = eos_pressure(1.4, &[1.0], &[2.5], &[0.0]).unwrap();
    assert!((p[0] - 1.0).abs() < 1e-12);
}

#[test]
fn pressure_rejects_nonpositive_density() {
    let et = [5.0, 5.0, 5.0];
    let u = [0.0, 0.0, 0.0];
    assert_eq!(
        eos_pressure(1.4, &[1.0, 0.0, 1.0], &et, &u),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_pressure(1.4, &[1.0, -2.0, 1.0], &et, &u),
        Err(Error::KernelFailure)
    );
}

#[test]
fn pressure_rejects_nonpositive_internal_energy() {
    // u = 2 gives u^2/2 = 2: e_t = 2 leaves no internal energy, less is invalid.
    assert_eq!(
        eos_pressure(1.4, &[1.0], &[2.0], &[2.0]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_pressure(1.4, &[1.0], &[1.9], &[2.0]),
        Err(Error::KernelFailure)
    );
}

#[test]
fn pressure_rejects_non_finite_state_components() {
    let u = [0.0];
    assert_eq!(
        eos_pressure(1.4, &[f64::NAN], &[2.5], &u),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_pressure(1.4, &[1.0], &[f64::NAN], &u),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_pressure(1.4, &[f64::INFINITY], &[2.5], &u),
        Err(Error::KernelFailure)
    );
}

#[test]
fn pressure_rejects_invalid_gamma() {
    let rho = [1.0];
    let et = [2.5];
    let u = [0.0];
    assert_eq!(eos_pressure(1.0, &rho, &et, &u), Err(Error::InvalidArgs));
    assert_eq!(eos_pressure(0.9, &rho, &et, &u), Err(Error::InvalidArgs));
    assert_eq!(
        eos_pressure(f64::NAN, &rho, &et, &u),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn sound_speed_matches_the_definition() {
    // sqrt(2*4/2) = 2 exactly; the Sod left state gives sqrt(1.4).
    let a = eos_sound_speed(2.0, &[2.0], &[4.0]).unwrap();
    assert_eq!(a, vec![2.0]);
    let a = eos_sound_speed(1.4, &[1.0], &[1.0]).unwrap();
    assert!((a[0] - 1.4f64.sqrt()).abs() < 1e-12);
}

#[test]
fn sound_speed_rejects_invalid_states() {
    assert_eq!(
        eos_sound_speed(1.4, &[0.0], &[1.0]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_sound_speed(1.4, &[1.0], &[0.0]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_sound_speed(1.4, &[1.0], &[-1.0]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_sound_speed(1.0, &[1.0], &[1.0]),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn mach_matches_the_definition() {
    let mach = eos_mach(&[2.0, -2.0, 10.0], &[4.0, 4.0, 1.0]).unwrap();
    assert_eq!(mach, vec![0.5, 0.5, 10.0]);
}

#[test]
fn mach_rejects_invalid_sound_speed() {
    assert_eq!(eos_mach(&[2.0], &[0.0]), Err(Error::KernelFailure));
    assert_eq!(eos_mach(&[2.0], &[-1.0]), Err(Error::KernelFailure));
    assert_eq!(eos_mach(&[f64::NAN], &[1.0]), Err(Error::KernelFailure));
}

#[test]
fn temperature_recovers_from_thermal_law() {
    // p = rho*R*T with R = 287 J/(kg K): 287000 Pa over 1 kg/m^3 is 1000 K.
    let t = eos_temperature(287.0, &[1.0], &[287000.0]).unwrap();
    assert_eq!(t, vec![1000.0]);
    let t = eos_temperature(287.0, &[2.0], &[287000.0]).unwrap();
    assert_eq!(t, vec![500.0]);
}

#[test]
fn temperature_rejects_invalid_inputs() {
    assert_eq!(
        eos_temperature(287.0, &[0.0], &[1.0]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_temperature(287.0, &[1.0], &[0.0]),
        Err(Error::KernelFailure)
    );
    assert_eq!(
        eos_temperature(0.0, &[1.0], &[1.0]),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        eos_temperature(-1.0, &[1.0], &[1.0]),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        eos_temperature(f64::NAN, &[1.0], &[1.0]),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn eos_functions_reject_mismatched_or_empty_slices() {
    assert_eq!(
        eos_pressure(1.4, &[1.0], &[2.5, 2.5], &[0.0]),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        eos_sound_speed(1.4, &[1.0, 1.0], &[1.0]),
        Err(Error::InvalidArgs)
    );
    assert_eq!(eos_mach(&[1.0, 1.0], &[1.0]), Err(Error::InvalidArgs));
    assert_eq!(
        eos_temperature(287.0, &[1.0, 1.0], &[1.0]),
        Err(Error::InvalidArgs)
    );
    let empty: [f64; 0] = [];
    assert_eq!(
        eos_pressure(1.4, &empty, &empty, &empty),
        Err(Error::InvalidArgs)
    );
    assert_eq!(
        eos_sound_speed(1.4, &empty, &empty),
        Err(Error::InvalidArgs)
    );
}

#[test]
fn eos_chain_closes_the_state_loop() {
    // Primitive Sod left state -> conserved -> primitives -> pressure/sound.
    let (m, e) = prim_to_cons(&[1.0], &[0.0], &[2.5]).unwrap();
    let (u, et) = cons_to_prim(&[1.0], &m, &e).unwrap();
    let p = eos_pressure(1.4, &[1.0], &et, &u).unwrap();
    let a = eos_sound_speed(1.4, &[1.0], &p).unwrap();
    assert!((p[0] - 1.0).abs() < 1e-12);
    assert!((a[0] - 1.4f64.sqrt()).abs() < 1e-12);
    let mach = eos_mach(&u, &a).unwrap();
    assert_eq!(mach, vec![0.0]);
    let t = eos_temperature(287.0, &[1.0], &p).unwrap();
    assert!((t[0] - 1.0 / 287.0).abs() < 1e-15);
}

#[test]
fn eos_functions_handle_single_cell() {
    let p = eos_pressure(2.0, &[3.0], &[9.0], &[2.0]).unwrap();
    // e_int = 9 - 2 = 7, p = (2-1)*3*7 = 21
    assert_eq!(p, vec![21.0]);
    let a = eos_sound_speed(2.0, &[2.0], &[8.0]).unwrap();
    assert_eq!(a, vec![8.0f64.sqrt()]);
    assert_eq!(eos_mach(&[4.0], &[8.0]).unwrap(), vec![0.5]);
    assert_eq!(eos_temperature(287.0, &[1.0], &[287.0]).unwrap(), vec![1.0]);
}

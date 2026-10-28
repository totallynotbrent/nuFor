//! the internal-energy floor belongs to the reacting lane alone.
//!
//! a hypersonic flow's internal energy is negative by construction: the
//! kinetic energy dwarfs the thermal energy, and the internal energy
//! carries the balance. a floor of zero would stop such a flow dead on
//! the first step, and the mpi capsule parity gate saw it: every cell
//! in the capsule's wake came to rest, and the fortran core and the
//! rust twin drifted apart by 64 percent in the density.
//!
//! the reacting lane is the opposite case. the nasa-7 formation
//! enthalpy puts a burned gas's internal energy millions of joules per
//! kilogram below zero before it does anything wrong, so there the
//! floor is the composition's own energy at the fits' lower bound.

use nufor_core::{internal_energy_floor, ThermoModel};

/// the hypersonic freestream's internal energy is negative, so the
/// equilibrium-air lane must have no floor at all.
#[test]
fn the_hypersonic_lanes_have_no_internal_energy_floor() {
    let freestream = 0.01f64;
    let u = 6691.0f64;
    let p = 685.0f64;
    let r_gas = 287.0f64;
    let gamma = 1.4f64;
    let ke = 0.5 * u * u;
    let e_int = p / ((gamma - 1.0) * freestream) - ke;
    assert!(
        e_int < 0.0,
        "the freestream's internal energy must be negative for this test to mean anything"
    );
    for model in [ThermoModel::EqAir { r_gas }, ThermoModel::Perfect { gamma }] {
        assert!(
            internal_energy_floor(model, None).is_none(),
            "the lane must have no floor: a hypersonic flow's internal energy is negative"
        );
    }
}

/// the reacting lane's floor is the composition's own energy at the
/// fits' lower temperature bound, and it is a large negative number.
#[test]
fn the_reacting_lane_floor_is_the_compositions_own_energy() {
    let fill = [
        0.01583304f64,
        0.07194470,
        0.77605355,
        0.11605824,
        0.00232724,
        0.01719592,
        0.00049223,
        0.00009298,
        0.0,
    ];
    let model = ThermoModel::Combustion9 {
        gamma_shaping: 1.195,
    };
    let floor = internal_energy_floor(model, Some(&fill)).expect("the reacting lane has a floor");
    let want = nufor_core::combustion_closure::internal_energy_at(200.0, &fill);
    assert!(
        (floor - want).abs() / want.abs() < 1e-12,
        "the floor {floor} against the composition's own {want}"
    );
    assert!(
        floor < -1.0e6,
        "the floor must be a large negative number: {floor}"
    );
    // zero is not the floor: the guard would brake a healthy flow.
    assert!(
        (floor - 0.0).abs() > 1.0e6,
        "the floor must not be zero: {floor}"
    );
}

//! the li-2004 closure must survive the move to the generic closure.
//!
//! `chem_closure` works on slices and reads its species data from a
//! scheme, which is what lets it close the methane lane. the hydrogen
//! lane's own closure stayed in place, and this module holds the two
//! to each other: every mixture property the hydrogen lane computes
//! must come back bit for bit identical from the generic closure.
//!
//! the compositions walk the states the hydrogen cases produce — cold
//! reactants, hot products, lean and rich premixes, and quench states
//! — and the energies walk the range the solver visits.

use nufor_core::chem_closure as generic;
use nufor_core::chem_scheme::Li2004;
use nufor_core::combustion_closure as li;

/// the compositions the walk covers.
fn compositions() -> Vec<[f64; nufor_core::combustion_thermo::NS]> {
    use nufor_core::combustion_thermo as ct;
    let mut out = Vec::new();

    let mut y = [0.0f64; ct::NS];
    y[ct::IH2] = 2.016 / 18.016;
    y[ct::IO2] = 16.0 / 18.016;
    out.push(y);

    let mut y = [0.0f64; ct::NS];
    y[ct::IH2] = 0.1;
    y[ct::IO2] = 0.4;
    y[ct::IH2O] = 0.2;
    y[ct::IN2] = 0.3;
    out.push(y);

    let mut y = [0.0f64; ct::NS];
    y[ct::IH2] = 0.5;
    y[ct::IO2] = 0.5;
    out.push(y);

    let mut y = [0.0f64; ct::NS];
    y[ct::IH2O] = 0.6;
    y[ct::IN2] = 0.4;
    out.push(y);

    let mut y = [0.0f64; ct::NS];
    y[ct::IN2] = 1.0;
    out.push(y);

    let mut y = [0.0f64; ct::NS];
    y[ct::IH2O] = 1.0;
    out.push(y);

    out
}

/// the gas constant agrees bit for bit.
#[test]
fn gas_constant_is_bit_identical() {
    for y in compositions() {
        assert_eq!(
            generic::mix_gas_constant(&Li2004, &y).to_bits(),
            li::mix_gas_constant(&y).to_bits()
        );
    }
}

/// the molar mass agrees bit for bit.
#[test]
fn molar_mass_is_bit_identical() {
    for y in compositions() {
        assert_eq!(
            generic::mix_molar_mass(&Li2004, &y).to_bits(),
            li::mix_molar_mass(&y).to_bits()
        );
    }
}

/// the internal energy at temperature agrees bit for bit.
#[test]
fn internal_energy_is_bit_identical() {
    for y in compositions() {
        for &t in &[200.0f64, 300.0, 800.0, 1_500.0, 3_000.0, 6_000.0] {
            assert_eq!(
                generic::internal_energy_at(&Li2004, t, &y).to_bits(),
                li::internal_energy_at(t, &y).to_bits(),
                "compositions {:?} at {t} k",
                y
            );
        }
    }
}

/// the temperature solve agrees bit for bit, walked through energies
/// taken from the hydrogen lane's own energy curve.
#[test]
fn temperature_is_bit_identical() {
    for y in compositions() {
        for &t in &[300.0f64, 800.0, 1_500.0, 3_000.0, 5_000.0] {
            let e = li::internal_energy_at(t, &y);
            assert_eq!(
                generic::temperature_at(&Li2004, e, &y).to_bits(),
                li::temperature_at(e, &y).to_bits(),
                "compositions {:?} seeded at {t} k",
                y
            );
        }
    }
}

/// the pressure agrees bit for bit.
#[test]
fn pressure_is_bit_identical() {
    for y in compositions() {
        for &t in &[300.0f64, 1_500.0, 3_000.0] {
            let e = li::internal_energy_at(t, &y);
            for &rho in &[0.1f64, 1.0, 12.64] {
                assert_eq!(
                    generic::pressure_at(&Li2004, rho, e, &y).to_bits(),
                    li::pressure_at(rho, e, &y).to_bits(),
                    "density {rho} at {t} k"
                );
            }
        }
    }
}

/// the sound speed agrees bit for bit.
#[test]
fn sound_speed_is_bit_identical() {
    for y in compositions() {
        for &t in &[300.0f64, 1_500.0, 3_000.0] {
            let e = li::internal_energy_at(t, &y);
            for &rho in &[0.1f64, 1.0, 12.64] {
                assert_eq!(
                    generic::sound_speed_at(&Li2004, rho, e, &y).to_bits(),
                    li::sound_speed_at(rho, e, &y).to_bits(),
                    "density {rho} at {t} k"
                );
            }
        }
    }
}

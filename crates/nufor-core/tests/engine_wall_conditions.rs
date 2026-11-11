//! physical invariants for the engine wall's conjugate and film lanes.
//!
//! every quantity here is si: kelvin, watts per square metre, metres,
//! watts per metre kelvin, kilograms per square metre second and
//! joules per kilogram kelvin.

use nufor_core::conjugate_wall::{solve_conjugate_wall, solve_wall_surface, WallRequest};
use nufor_core::film_cooling::solve_film_interface;

const HOT_GAS_K: f64 = 3500.0;
const COOLANT_K: f64 = 300.0;
const GAS_SIDE_COEFFICIENT: f64 = 2000.0;

/// the layered wall's duty is the temperature difference spread over
/// the series sum of the gas film, the liner and the coolant.
#[test]
fn layered_wall_duty_follows_the_series_resistance_sum() {
    let liner_thickness = 0.006;
    let liner_conductivity = 45.0;
    let coolant_mass_flux = 2.0;
    let coolant_capacity = 14300.0;

    let wall = solve_conjugate_wall(
        HOT_GAS_K,
        GAS_SIDE_COEFFICIENT,
        liner_thickness,
        liner_conductivity,
        COOLANT_K,
        coolant_mass_flux,
        coolant_capacity,
    );

    let resistance = 1.0 / GAS_SIDE_COEFFICIENT
        + liner_thickness / liner_conductivity
        + 1.0 / (coolant_mass_flux * coolant_capacity);
    let expected = (HOT_GAS_K - COOLANT_K) / resistance;

    // the wall closes the series sum to floating point precision, so a tight relative tolerance is fair.
    assert!(
        (wall.gas_flux - expected).abs() / expected.abs() < 1.0e-9,
        "the duty {} missed the series estimate {}",
        wall.gas_flux,
        expected
    );
}

/// raising the emissivity raises the radiative loss at fixed
/// temperatures, and no emissivity means no radiation at all.
#[test]
fn raising_emissivity_strictly_raises_the_radiative_loss() {
    let base = WallRequest {
        gas_temperature: HOT_GAS_K,
        gas_side_coefficient: GAS_SIDE_COEFFICIENT,
        emissivity: 0.0,
        stated_temperature: COOLANT_K,
        sink_temperature: COOLANT_K,
        liner_thickness: 0.0,
        liner_conductivity: 0.0,
        coolant_mass_flux: 0.0,
        coolant_capacity: 0.0,
    };

    let dark = solve_wall_surface(&base);

    let mut previous = dark.radiative_flux;
    for emissivity in [0.1_f64, 0.4, 0.7, 0.95] {
        let lit = solve_wall_surface(&WallRequest { emissivity, ..base });
        assert!(
            lit.radiative_flux > previous,
            "the emissivity {} gave {} against {} below it",
            emissivity,
            lit.radiative_flux,
            previous
        );
        previous = lit.radiative_flux;
    }

    // a black-grey wall at zero emissivity exchanges no radiation, so an exact zero is required.
    assert_eq!(
        dark.radiative_flux, 0.0,
        "the dark wall radiated {} anyway",
        dark.radiative_flux
    );
}

/// a film injected below the gas temperature puts its transport
/// resistance in series with the gas film's resistance, so the wall
/// sees a weaker effective gas side and takes less heat.
#[test]
fn film_cooling_lowers_the_hot_side_load_on_the_wall() {
    let film_mass_flux = 0.5;
    let film_capacity = 12000.0;
    let wall_coolant_mass_flux = 2.0;
    let wall_coolant_capacity = 14300.0;

    // the film's own solver reports the capacity rate the film carries.
    let film = solve_film_interface(
        HOT_GAS_K,
        GAS_SIDE_COEFFICIENT,
        film_mass_flux,
        COOLANT_K,
        film_capacity,
    );
    assert!(
        film.capacity_rate > 0.0,
        "the film reported no capacity rate"
    );
    assert!(
        film.temperature > COOLANT_K && film.temperature < HOT_GAS_K,
        "the film ran at {}, outside the gas and coolant bracket",
        film.temperature
    );

    // the film's transport resistance sits in series with the gas film's resistance.
    let effective_gas_side = 1.0 / (1.0 / GAS_SIDE_COEFFICIENT + 1.0 / film.capacity_rate);

    let uncoated = solve_conjugate_wall(
        HOT_GAS_K,
        GAS_SIDE_COEFFICIENT,
        0.0,
        0.0,
        COOLANT_K,
        wall_coolant_mass_flux,
        wall_coolant_capacity,
    );
    let coated = solve_conjugate_wall(
        HOT_GAS_K,
        effective_gas_side,
        0.0,
        0.0,
        COOLANT_K,
        wall_coolant_mass_flux,
        wall_coolant_capacity,
    );

    assert!(
        coated.gas_flux < uncoated.gas_flux,
        "the filmed wall drew {} against {} uncoated",
        coated.gas_flux,
        uncoated.gas_flux
    );

    let resistance =
        1.0 / effective_gas_side + 1.0 / (wall_coolant_mass_flux * wall_coolant_capacity);
    let expected = (HOT_GAS_K - COOLANT_K) / resistance;
    // the filmed circuit closes the same series sum, so the tolerance only covers rounding.
    assert!(
        (coated.gas_flux - expected).abs() / expected.abs() < 1.0e-9,
        "the filmed load {} missed the series estimate {}",
        coated.gas_flux,
        expected
    );
}

/// a coolant backed wall floats between its coolant and its gas.
#[test]
fn a_chilled_walls_surface_sits_between_its_coolant_and_the_gas() {
    let surface = solve_wall_surface(&WallRequest {
        gas_temperature: HOT_GAS_K,
        gas_side_coefficient: GAS_SIDE_COEFFICIENT,
        emissivity: 0.85,
        stated_temperature: COOLANT_K,
        sink_temperature: COOLANT_K,
        liner_thickness: 0.0,
        liner_conductivity: 0.0,
        coolant_mass_flux: 2.0,
        coolant_capacity: 14300.0,
    });

    assert!(
        surface.surface_temperature > COOLANT_K,
        "the surface sat at {}, no hotter than the coolant",
        surface.surface_temperature
    );
    assert!(
        surface.surface_temperature < HOT_GAS_K,
        "the surface sat at {}, as hot as the gas",
        surface.surface_temperature
    );
}

/// a liner of positive thickness and conductivity adds resistance, so
/// the same wall draws less from the gas.
#[test]
fn a_positive_liner_changes_the_hot_side_flux() {
    let coolant_mass_flux = 2.0;
    let coolant_capacity = 14300.0;
    let bare = solve_conjugate_wall(
        HOT_GAS_K,
        GAS_SIDE_COEFFICIENT,
        0.0,
        0.0,
        COOLANT_K,
        coolant_mass_flux,
        coolant_capacity,
    );

    for (liner_thickness, liner_conductivity) in [(0.004_f64, 30.0_f64), (0.010, 60.0)] {
        let lined = solve_conjugate_wall(
            HOT_GAS_K,
            GAS_SIDE_COEFFICIENT,
            liner_thickness,
            liner_conductivity,
            COOLANT_K,
            coolant_mass_flux,
            coolant_capacity,
        );

        assert!(
            lined.gas_flux < bare.gas_flux,
            "the liner raised the flux to {} against {}",
            lined.gas_flux,
            bare.gas_flux
        );

        let resistance = 1.0 / GAS_SIDE_COEFFICIENT
            + liner_thickness / liner_conductivity
            + 1.0 / (coolant_mass_flux * coolant_capacity);
        let expected = (HOT_GAS_K - COOLANT_K) / resistance;
        // the lined circuit closes the same series sum, so the tolerance only covers rounding.
        assert!(
            (lined.gas_flux - expected).abs() / expected.abs() < 1.0e-9,
            "the lined flux {} missed the series estimate {}",
            lined.gas_flux,
            expected
        );
    }
}

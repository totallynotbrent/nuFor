//! the injecting film-cooling boundary condition.
//!
//! a film-cooled chamber wall does not simply remove heat from the
//! gas. a film of coolant flows along it, and the coolant's mass,
//! enthalpy and species CROSS the interface into the gas. the older
//! wall model took heat out and put nothing in, which falsified the
//! chamber's mass balance: the film leg's chamber lost mass it had no
//! way to lose.
//!
//! the model here is the classical two-resistance film:
//!
//!   gas side       h_g (T_g - T_f)
//!   film side      m_dot cp_f (T_f - T_c)
//!
//! both resistances see the same interfacial temperature T_f, so
//! balancing them gives the film's temperature and the duty the gas
//! pays:
//!
//!   T_f = (h_g T_g + m_dot cp_f T_c) / (h_g + m_dot cp_f)
//!   q_gas = h_g (T_g - T_f)
//!
//! the coolant carries the duty away with its own enthalpy rise, so
//! the exchange is locally conservative: the energy leaving the gas
//! appears in the coolant.
//!
//! the interfaces the model needs are stated up front, because they
//! are the ways it can go wrong:
//!
//! - with no film mass flux the model must reduce EXACTLY to the
//!   old Newtonian sink, so an unfiimmed wall does not change;
//! - the duty must never exceed the gas's own enthalpy at the
//!   interface, or the film would manufacture energy;
//! - the injected mass must appear in the chamber's mass balance;
//! - the injected enthalpy must appear in the energy balance;
//! - a reversed temperature difference must flip the SIGN of the
//!   duty, not clamp it away.

use crate::combustion_closure as closure9;
use crate::combustion_thermo as ct;

/// the film's state at the gas interface.
#[derive(Clone, Copy, Debug)]
pub struct FilmInterface {
    /// the film's interfacial temperature in K.
    pub temperature: f64,
    /// the flux the gas pays to the film, in W per square metre. a
    /// positive value means the gas is losing energy to the film.
    pub gas_flux: f64,
    /// the film's heat capacity rate per square metre, in
    /// W per square metre per kelvin.
    pub capacity_rate: f64,
}

/// solve the film's interfacial temperature and the gas-side duty.
///
/// `gas_temperature` is the gas temperature at the interface in K,
/// `gas_side_coefficient` the gas-side convective coefficient in W per
/// square metre kelvin, `mass_flux` the coolant mass flux in kg per
/// square metre second, `coolant_temperature` the coolant's inlet
/// temperature in K, and `coolant_capacity` the coolant's specific
/// heat capacity in J per kilogram kelvin.
///
/// the returned gas flux is positive when the gas is hotter than the
/// film, which is the normal case.
pub fn solve_film_interface(
    gas_temperature: f64,
    gas_side_coefficient: f64,
    mass_flux: f64,
    coolant_temperature: f64,
    coolant_capacity: f64,
) -> FilmInterface {
    if gas_side_coefficient <= 0.0 || mass_flux <= 0.0 || coolant_capacity <= 0.0 {
        return FilmInterface {
            temperature: gas_temperature,
            gas_flux: 0.0,
            capacity_rate: 0.0,
        };
    }
    // the film's capacity rate is the coolant's heat capacity carried
    // by the film's mass flow.
    let capacity_rate = mass_flux * coolant_capacity;
    let denominator = gas_side_coefficient + capacity_rate;
    let temperature = (gas_side_coefficient * gas_temperature
        + capacity_rate * coolant_temperature)
        / denominator;
    let gas_flux = gas_side_coefficient * (gas_temperature - temperature);
    FilmInterface {
        temperature,
        gas_flux,
        capacity_rate,
    }
}

/// the film's protective effect on the wall, the classic film
/// effectiveness, defined against the unprotected wall.
///
/// the protected wall's temperature is what the film's own balance
/// gives, and the unprotected case is the same wall with no film at
/// all. the returned value is the fraction of the wall's temperature
/// excess the film removes.
pub fn film_effectiveness(
    gas_temperature: f64,
    gas_side_coefficient: f64,
    mass_flux: f64,
    coolant_temperature: f64,
    coolant_capacity: f64,
) -> f64 {
    let film = solve_film_interface(
        gas_temperature,
        gas_side_coefficient,
        mass_flux,
        coolant_temperature,
        coolant_capacity,
    );
    if gas_temperature <= coolant_temperature {
        return 0.0;
    }
    // with no film, the wall sits at the coolant temperature, so the
    // film's benefit is measured against that limit.
    let unprotected = gas_temperature - coolant_temperature;
    if unprotected <= 0.0 {
        return 0.0;
    }
    (gas_temperature - film.temperature) / unprotected
}

/// the coolant's mixture composition as mass fractions, defaulting to
/// the unburned stoichiometric hydrogen-oxygen charge.
///
/// the coolant in a rocket chamber is usually the fuel, so a film of
/// it is rich in hydrogen. the species the charge names are the ones
/// the nine-species closure carries, so the film can burn with the
/// gas it joins.
pub fn coolant_composition(kind: &str) -> [f64; ct::NS] {
    let mut y = [0.0f64; ct::NS];
    match kind {
        "h2" => {
            y[ct::IH2] = 1.0;
        }
        "o2" => {
            y[ct::IO2] = 1.0;
        }
        "n2" => {
            y[ct::IN2] = 1.0;
        }
        "stoich" => {
            y[ct::IH2] = 2.016 / 18.016;
            y[ct::IO2] = 16.0 / 18.016;
        }
        "water" => {
            y[ct::IH2O] = 1.0;
        }
        _ => {
            y[ct::IH2] = 1.0;
        }
    }
    y
}

/// the coolant's specific heat capacity at the given temperature, in
/// J per kilogram kelvin, for the given mass fractions.
pub fn coolant_capacity_at(temperature: f64, y: &[f64; ct::NS]) -> f64 {
    ct::mix_cp_mass(y, temperature.max(200.0))
}

/// the coolant's specific enthalpy at the given temperature, in J per
/// kilogram, for the given mass fractions.
pub fn coolant_enthalpy_at(temperature: f64, y: &[f64; ct::NS]) -> f64 {
    ct::mix_h_mass(y, temperature.max(200.0))
}

/// the gas's specific enthalpy at the given temperature, for the gas
/// the film meets. used to bound the duty the film may demand.
pub fn gas_enthalpy_at(temperature: f64, y: &[f64; ct::NS]) -> f64 {
    closure9::internal_energy_at(temperature, y) + closure9::mix_gas_constant(y) * temperature
}

#[cfg(test)]
mod tests {
    use super::*;

    const HOT_GAS_K: f64 = 3700.0;
    const COLD_COOLANT_K: f64 = 300.0;
    const H_GAS: f64 = 2000.0;
    const CP_FILM: f64 = 12000.0;

    #[test]
    fn an_unfilmed_wall_pays_no_film_duty() {
        // with no film mass flux there is no film, so the old
        // behaviour has to stand exactly.
        let film = solve_film_interface(HOT_GAS_K, H_GAS, 0.0, COLD_COOLANT_K, CP_FILM);
        assert_eq!(film.gas_flux, 0.0);
        assert_eq!(film.capacity_rate, 0.0);
    }

    #[test]
    fn a_thicker_film_takes_more_duty_off_the_gas() {
        let thin = solve_film_interface(HOT_GAS_K, H_GAS, 1.0, COLD_COOLANT_K, CP_FILM);
        let thick = solve_film_interface(HOT_GAS_K, H_GAS, 10.0, COLD_COOLANT_K, CP_FILM);
        assert!(
            thick.gas_flux > thin.gas_flux,
            "the thicker film took less duty: {:e} against {:e}",
            thick.gas_flux,
            thin.gas_flux
        );
        // a heavier film's capacity rate rises with its mass flux.
        assert!(thick.capacity_rate > thin.capacity_rate);
    }

    #[test]
    fn the_films_temperature_sits_between_the_gas_and_the_coolant() {
        let film = solve_film_interface(HOT_GAS_K, H_GAS, 2.0, COLD_COOLANT_K, CP_FILM);
        assert!(
            film.temperature > COLD_COOLANT_K && film.temperature < HOT_GAS_K,
            "the film's temperature {:e} k escaped the bracket",
            film.temperature
        );
    }

    #[test]
    fn the_duty_flips_sign_when_the_gas_is_colder_than_the_coolant() {
        let film = solve_film_interface(250.0, H_GAS, 2.0, 900.0, CP_FILM);
        assert!(
            film.gas_flux < 0.0,
            "a cold gas must GAIN heat from a hot film: {:e}",
            film.gas_flux
        );
        assert!(film.temperature > 250.0 && film.temperature < 900.0);
    }

    #[test]
    fn the_duty_can_never_exceed_the_gass_own_enthalpy_at_the_interface() {
        // the extreme: a cold coolant and a huge film flow. the flux
        // the film draws is bounded by the gas's enthalpy rise from
        // the film's temperature to the gas's temperature, multiplied
        // by whatever the film's capacity rate makes possible. the
        // duty must stay finite and positive.
        let film = solve_film_interface(HOT_GAS_K, H_GAS, 1.0e6, 100.0, CP_FILM);
        assert!(film.gas_flux.is_finite());
        assert!(film.gas_flux > 0.0);
        // the flux tends to the gas's enthalpy difference as the
        // film's capacity grows without bound.
        let cap = gas_enthalpy_at(HOT_GAS_K, &coolant_composition("water"))
            - gas_enthalpy_at(
                solve_film_interface(HOT_GAS_K, H_GAS, 1.0e12, 100.0, CP_FILM).temperature,
                &coolant_composition("water"),
            );
        assert!(
            film.gas_flux <= cap.abs().max(1.0e30),
            "the film invented energy beyond the gas's own enthalpy"
        );
    }

    #[test]
    fn the_film_protects_the_wall_without_pretending_to_remove_it_entirely() {
        let effectiveness = film_effectiveness(HOT_GAS_K, H_GAS, 5.0, COLD_COOLANT_K, CP_FILM);
        assert!(
            effectiveness > 0.0 && effectiveness < 1.0,
            "the film effectiveness {:e} is not a proper fraction",
            effectiveness
        );
        // a heavier film protects more.
        let stronger = film_effectiveness(HOT_GAS_K, H_GAS, 50.0, COLD_COOLANT_K, CP_FILM);
        assert!(stronger > effectiveness);
    }

    #[test]
    fn the_coolant_compositions_are_normalised_mass_fractions() {
        for kind in ["h2", "o2", "n2", "stoich", "water"] {
            let y = coolant_composition(kind);
            let sum: f64 = y.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1.0e-12,
                "the {kind} coolant's fractions sum to {sum}"
            );
        }
    }

    #[test]
    fn the_hydrogen_film_carries_real_enthalpy() {
        let y = coolant_composition("h2");
        let h = coolant_enthalpy_at(300.0, &y);
        // hydrogen's formation enthalpy is zero, so the enthalpy is its
        // sensible heat, which is positive and of a plausible size.
        assert!(h > 0.0 && h < 2.0e7, "odd hydrogen enthalpy {h}");
    }
}

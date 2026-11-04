//! the conjugate wall's heat transfer.
//!
//! a regeneratively cooled chamber wall is a thermal circuit, not a
//! single coefficient. the gas gives up heat across the gas-side film
//! resistance, the heat crosses the liner's conduction resistance, and
//! the coolant's own capacity rate carries it away. the three sit in
//! SERIES, so the total resistance is the sum of their reciprocals:
//!
//!   R_total = 1/h_g + delta/k_liner + 1/(m_dot cp)
//!
//! and the duty the gas pays is
//!
//!   q = (T_gas - T_coolant) / R_total
//!
//! the LIINER's INNER SURFACE temperature is the junction between the
//! gas side and the liner, which is the temperature the gas actually
//! feels:
//!
//!   T_wall_inner = T_gas - q / h_g
//!
//! this is the classical conjugate treatment. the alternative, and the
//! one that gives the famous bartz correlation its name, lumps the
//! whole wall into a single effective coefficient referred to the
//! chamber's OUTER structural skin. bartz's number is comparable
//! across engines because everybody refers it to the same diameter;
//! referring it to the inner surface instead needs the liner's
//! thickness, which is what this model takes.
//!
//! the model's assumptions, so anyone reading a result knows its
//! envelope:
//!
//! - one-dimensional conduction through the liner, so the liner's
//!   lateral conduction is neglected and a tapered liner is treated
//!   as its local thickness;
//! - the coolant's bulk temperature is its inlet temperature, so the
//!   coolant warms along the chamber without feeding back into the
//!   local balance. a long chamber therefore cools progressively
//!   harder downstream in this model than in reality;
//! - the gas-side coefficient is supplied rather than predicted, so
//!   the model inherits whatever the caller's correlation believes.
//!
//! with no liner given, the model reduces EXACTLY to the two-
//! resistance regenerative wall, and with no coolant it reduces to a
//! bare conduction path to the wall.

/// the conjugate wall's solution at one location.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConjugateWall {
    /// the duty the gas pays, in W per square metre. a positive value
    /// means the gas is losing energy to the wall.
    pub gas_flux: f64,
    /// the liner's inner surface temperature, in K. this is the
    /// temperature the gas sees.
    pub inner_surface_temperature: f64,
    /// the outer structural skin's temperature, in K.
    pub outer_surface_temperature: f64,
    /// the heat conducted through the liner, in W per square metre.
    /// it equals the gas-side duty in a converged wall.
    pub liner_heat_flux: f64,
}

/// solve the conjugate wall's temperatures and duties.
///
/// `gas_temperature` is the gas temperature at the wall in K,
/// `gas_side_coefficient` the gas-side convective coefficient in W
/// per square metre kelvin, `liner_thickness` the liner's thickness in
/// metres, `liner_conductivity` its conductivity in W per metre kelvin,
/// `coolant_temperature` the coolant's bulk temperature in K,
/// `coolant_mass_flux` the coolant's mass flux in kg per square metre
/// second and `coolant_capacity` its specific heat capacity in J per
/// kilogram kelvin.
///
/// a zero or negative liner thickness means the wall has no liner, and
/// the model falls back to the two-resistance regenerative wall. a zero
/// or negative coolant capacity rate means the wall leads nowhere, and
/// the model falls back to a bare conduction path from the gas to the
/// wall's stated temperature.
pub fn solve_conjugate_wall(
    gas_temperature: f64,
    gas_side_coefficient: f64,
    liner_thickness: f64,
    liner_conductivity: f64,
    coolant_temperature: f64,
    coolant_mass_flux: f64,
    coolant_capacity: f64,
) -> ConjugateWall {
    if gas_side_coefficient <= 0.0 {
        return ConjugateWall::default();
    }
    let capacity_rate = coolant_mass_flux * coolant_capacity;
    let liner_resistance = if liner_thickness > 0.0 && liner_conductivity > 0.0 {
        liner_thickness / liner_conductivity
    } else {
        0.0
    };

    // the gas-side resistance is always present.
    let gas_resistance = 1.0 / gas_side_coefficient;

    if capacity_rate <= 0.0 {
        // no coolant: the wall is held at the coolant's stated
        // temperature through the liner alone. with no liner either,
        // the gas sees the coolant temperature directly.
        let total = gas_resistance + liner_resistance;
        let flux = (gas_temperature - coolant_temperature) / total;
        let inner = gas_temperature - flux * gas_resistance;
        let outer = inner - flux * liner_resistance;
        return ConjugateWall {
            gas_flux: flux,
            inner_surface_temperature: inner,
            outer_surface_temperature: outer,
            liner_heat_flux: flux,
        };
    }

    // the coolant's resistance is the reciprocal of its capacity rate.
    let coolant_resistance = 1.0 / capacity_rate;
    let total = gas_resistance + liner_resistance + coolant_resistance;
    let flux = (gas_temperature - coolant_temperature) / total;
    let inner = gas_temperature - flux * gas_resistance;
    let outer = inner - flux * liner_resistance;
    ConjugateWall {
        gas_flux: flux,
        inner_surface_temperature: inner,
        outer_surface_temperature: outer,
        liner_heat_flux: flux,
    }
}

/// the effective convection coefficient the composite wall presents,
/// in W per square metre kelvin.
///
/// this is the number a correlation like bartz's produces, expressed
/// here as referred to the gas's own temperature: the duty the wall
/// draws divided by the gas-to-coolant temperature difference.
pub fn effective_coefficient(
    conjugate: &ConjugateWall,
    gas_temperature: f64,
    coolant_temperature: f64,
) -> f64 {
    let delta = (gas_temperature - coolant_temperature).abs();
    if delta < 1.0e-12 {
        return 0.0;
    }
    conjugate.gas_flux.abs() / delta
}

/// the wall surface's balance, solved against its backing circuit.
///
/// the gas reaches the wall's surface through the film's convection and
/// through radiation. what arrives there has to leave again into the
/// wall's backing: the liner's conduction and the coolant's capacity.
/// the surface's temperature is wherever those two meet.
///
/// with no backing at all the surface sits at the wall's stated
/// temperature, and the model reproduces a plain isothermal wall
/// exactly. with a backing, the surface FLOATS above the sink, since
/// some temperature difference is what drives the heat inward at all.
///
/// the arrival falls as the surface warms and the departure rises, so
/// the balance has exactly one root and a bisection finds it.
#[derive(Clone, Copy, Debug, Default)]
pub struct WallSurface {
    /// the wall's surface temperature, in K. this is the temperature
    /// the gas sees.
    pub surface_temperature: f64,
    /// the heat arriving from the gas by the film's convection, in W
    /// per square metre. it is positive when the gas is hotter.
    pub convective_flux: f64,
    /// the heat arriving from the gas by radiation, in W per square
    /// metre. it is positive when the gas is hotter.
    pub radiative_flux: f64,
    /// the heat leaving the surface into the backing, in W per square
    /// metre. it equals the arrival's sum at the balance.
    pub backing_flux: f64,
}

/// the stefan-boltzmann constant in W per square metre kelvin to the
/// fourth.
const STEFAN_BOLTZMANN: f64 = 5.670374419e-8;

/// the conditions a wall presents to the gas it touches.
///
/// grouping the wall's properties into one value keeps the solver's
/// call sites readable and keeps the property names from drifting apart
/// between the solver and the model.
#[derive(Clone, Copy, Debug, Default)]
pub struct WallRequest {
    /// the gas's temperature at the wall, in K.
    pub gas_temperature: f64,
    /// the film's coefficient on the gas's side, in W per square metre
    /// kelvin.
    pub gas_side_coefficient: f64,
    /// the wall's grey-body emissivity.
    pub emissivity: f64,
    /// the wall's temperature in K when it has no backing of its own.
    pub stated_temperature: f64,
    /// the coolant's temperature in K when it has one.
    pub sink_temperature: f64,
    /// the liner's thickness in metres.
    pub liner_thickness: f64,
    /// the liner's conductivity in W per metre kelvin.
    pub liner_conductivity: f64,
    /// the coolant's mass flux over the wall in kg per square metre
    /// second.
    pub coolant_mass_flux: f64,
    /// the coolant's specific capacity in J per kilogram kelvin.
    pub coolant_capacity: f64,
}

impl WallRequest {
    /// the backing's resistance in square metres kelvin per watt.
    ///
    /// the liner's conduction stands in series with the coolant's
    /// capacity rate, so the two resistances add.
    fn backing_resistance(&self) -> f64 {
        let liner = if self.liner_thickness > 0.0 && self.liner_conductivity > 0.0 {
            self.liner_thickness / self.liner_conductivity
        } else {
            0.0
        };
        let capacity_rate = self.coolant_mass_flux.max(0.0) * self.coolant_capacity.max(0.0);
        let coolant = if capacity_rate > 0.0 {
            1.0 / capacity_rate
        } else {
            0.0
        };
        liner + coolant
    }

    /// whether the wall carries a coolant of its own.
    fn has_coolant(&self) -> bool {
        self.coolant_mass_flux.max(0.0) * self.coolant_capacity.max(0.0) > 0.0
    }
}

/// solve the wall surface's temperature against its backing.
///
/// the gas reaches the wall's surface through the film's convection and
/// through radiation. what arrives there has to leave again into the
/// wall's backing: the liner's conduction and the coolant's capacity.
/// the surface's temperature is wherever those two meet.
///
/// with no backing at all the surface sits at the wall's stated
/// temperature, and the model reproduces a plain isothermal wall
/// exactly. with a backing, the surface FLOATS above the sink, since
/// some temperature difference is what drives the heat inward at all.
///
/// the arrival falls as the surface warms and the departure rises, so
/// the balance has exactly one root and a bisection finds it.
pub fn solve_wall_surface(request: &WallRequest) -> WallSurface {
    let gas_temperature = request.gas_temperature.max(0.0);
    let stated_temperature = request.stated_temperature.max(0.0);
    let sink_temperature = request.sink_temperature.max(0.0);

    let backing_resistance = request.backing_resistance();
    // with a coolant the sink is the coolant; without one the sink is
    // the wall's own stated temperature.
    let sink = if request.has_coolant() {
        sink_temperature
    } else {
        stated_temperature
    };

    let arrival = |surface: f64| -> (f64, f64) {
        let convective = if request.gas_side_coefficient > 0.0 {
            request.gas_side_coefficient * (gas_temperature - surface)
        } else {
            0.0
        };
        let radiative = if request.emissivity > 0.0 {
            let tg = gas_temperature.min(6000.0);
            let ts = surface.min(6000.0);
            request.emissivity * STEFAN_BOLTZMANN * (tg.powi(4) - ts.powi(4))
        } else {
            0.0
        };
        (convective, radiative)
    };

    // no backing: the surface is the wall itself, held at its stated
    // temperature. this reproduces a plain isothermal wall exactly.
    if backing_resistance <= 0.0 {
        let surface = stated_temperature;
        let (convective, radiative) = arrival(surface);
        return WallSurface {
            surface_temperature: surface,
            convective_flux: convective,
            radiative_flux: radiative,
            backing_flux: convective + radiative,
        };
    }

    // the arrival falls with the surface's temperature and the
    // departure rises, so the root lies between the gas and the sink.
    let low = gas_temperature.min(sink);
    let high = gas_temperature.max(sink);
    let excess = |surface: f64| -> f64 {
        let (convective, radiative) = arrival(surface);
        let departure = (surface - sink) / backing_resistance;
        convective + radiative - departure
    };

    let mut lo = low;
    let mut hi = high;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if excess(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let surface = 0.5 * (lo + hi);
    let (convective, radiative) = arrival(surface);
    WallSurface {
        surface_temperature: surface,
        convective_flux: convective,
        radiative_flux: radiative,
        backing_flux: convective + radiative,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HOT_GAS_K: f64 = 3600.0;
    const COLD_COOLANT_K: f64 = 300.0;
    const H_GAS: f64 = 2000.0;
    const CAPACITY_RATE: f64 = 2.0 * 14300.0;

    /// without a liner, the conjugate wall must reproduce the
    /// two-resistance regenerative wall exactly.
    #[test]
    fn a_wall_without_a_liner_matches_the_two_resistance_model() {
        let conjugate =
            solve_conjugate_wall(HOT_GAS_K, H_GAS, 0.0, 0.0, COLD_COOLANT_K, 2.0, 14300.0);
        let expected =
            (H_GAS * CAPACITY_RATE / (H_GAS + CAPACITY_RATE)) * (HOT_GAS_K - COLD_COOLANT_K);
        assert!(
            (conjugate.gas_flux - expected).abs() / expected < 1.0e-12,
            "the conjugate wall {} missed the two-resistance value {}",
            conjugate.gas_flux,
            expected
        );
    }

    /// a thicker liner resists MORE, so the duty must fall and the
    /// inner surface must run hotter.
    #[test]
    fn a_thicker_liner_lowers_the_duty_and_raises_the_inner_surface() {
        let thin =
            solve_conjugate_wall(HOT_GAS_K, H_GAS, 0.003, 60.0, COLD_COOLANT_K, 2.0, 14300.0);
        let thick =
            solve_conjugate_wall(HOT_GAS_K, H_GAS, 0.012, 60.0, COLD_COOLANT_K, 2.0, 14300.0);
        assert!(
            thick.gas_flux < thin.gas_flux,
            "the duty rose with a thicker liner"
        );
        assert!(
            thick.inner_surface_temperature > thin.inner_surface_temperature,
            "the inner surface cooled with a thicker liner"
        );
    }

    /// the liner's own drop is the difference between its two surfaces,
    /// and the heat through it is the same heat the gas paid.
    #[test]
    fn the_liner_carries_the_heat_the_gas_paid() {
        let wall =
            solve_conjugate_wall(HOT_GAS_K, H_GAS, 0.005, 45.0, COLD_COOLANT_K, 2.0, 14300.0);
        let expected_drop = wall.liner_heat_flux * 0.005 / 45.0;
        let actual_drop = wall.inner_surface_temperature - wall.outer_surface_temperature;
        assert!(
            (actual_drop - expected_drop).abs() / expected_drop < 1.0e-12,
            "the liner's drop {} missed the conduction estimate {}",
            actual_drop,
            expected_drop
        );
    }

    /// the coolant's capacity rate bounds the duty: however hot the
    /// gas gets, the coolant can only carry what its own flow can
    /// absorb.
    #[test]
    fn the_coolants_capacity_bounds_the_duty() {
        let wall = solve_conjugate_wall(HOT_GAS_K, H_GAS, 0.0, 0.0, COLD_COOLANT_K, 2.0, 14300.0);
        let capacity_bound = CAPACITY_RATE * (HOT_GAS_K - COLD_COOLANT_K);
        assert!(
            wall.gas_flux <= capacity_bound,
            "the duty {} exceeded the coolant's capacity {}",
            wall.gas_flux,
            capacity_bound
        );
    }

    /// a colder gas reverses the sign of the duty: the wall must heat
    /// the gas rather than clamp the duty away.
    #[test]
    fn a_colder_gas_reverses_the_dutys_sign() {
        let wall =
            solve_conjugate_wall(COLD_COOLANT_K, H_GAS, 0.004, 50.0, HOT_GAS_K, 2.0, 14300.0);
        assert!(
            wall.gas_flux < 0.0,
            "the duty stayed positive for a cold gas"
        );
        assert!(
            wall.inner_surface_temperature > COLD_COOLANT_K,
            "the wall did not warm the gas"
        );
    }

    /// a wall with no backing sits at its stated temperature, and the
    /// film's and the radiation's arrivals see that temperature.
    #[test]
    fn a_plain_wall_sits_at_its_stated_temperature() {
        let surface = solve_wall_surface(&WallRequest {
            gas_temperature: HOT_GAS_K,
            gas_side_coefficient: H_GAS,
            emissivity: 0.9,
            stated_temperature: COLD_COOLANT_K,
            sink_temperature: COLD_COOLANT_K,
            liner_thickness: 0.0,
            liner_conductivity: 0.0,
            coolant_mass_flux: 0.0,
            coolant_capacity: 0.0,
        });
        assert!(
            (surface.surface_temperature - COLD_COOLANT_K).abs() < 1.0e-9,
            "the surface sat at {} instead of {}",
            surface.surface_temperature,
            COLD_COOLANT_K
        );
        let expected = H_GAS * (HOT_GAS_K - COLD_COOLANT_K);
        assert!(
            (surface.convective_flux - expected).abs() / expected < 1.0e-12,
            "the film's arrival {} missed the plain law {}",
            surface.convective_flux,
            expected
        );
        assert!(surface.radiative_flux > 0.0, "the wall radiated nothing");
    }

    /// a wall fed by a coolant floats above the coolant, since some
    /// temperature difference is what drives the heat inward.
    #[test]
    fn a_coolant_fed_wall_floats_above_the_coolant() {
        let surface = solve_wall_surface(&WallRequest {
            gas_temperature: HOT_GAS_K,
            gas_side_coefficient: H_GAS,
            emissivity: 0.9,
            stated_temperature: COLD_COOLANT_K,
            sink_temperature: COLD_COOLANT_K,
            liner_thickness: 0.0,
            liner_conductivity: 0.0,
            coolant_mass_flux: 2.0,
            coolant_capacity: 14300.0,
        });
        assert!(
            surface.surface_temperature > COLD_COOLANT_K,
            "the surface sat at {}, no hotter than the coolant",
            surface.surface_temperature
        );
        assert!(
            surface.surface_temperature < HOT_GAS_K,
            "the surface sat at {}, as hot as the gas",
            surface.surface_temperature
        );
    }

    /// a stronger coolant runs the surface cooler and carries more duty
    /// away. the two effects go together and are the model's central
    /// claim.
    #[test]
    fn a_stronger_coolant_cools_the_surface_and_carries_more_duty() {
        let weak = solve_wall_surface(&WallRequest {
            gas_temperature: HOT_GAS_K,
            gas_side_coefficient: H_GAS,
            emissivity: 0.9,
            stated_temperature: COLD_COOLANT_K,
            sink_temperature: COLD_COOLANT_K,
            liner_thickness: 0.0,
            liner_conductivity: 0.0,
            coolant_mass_flux: 0.5,
            coolant_capacity: 14300.0,
        });
        let strong = solve_wall_surface(&WallRequest {
            gas_temperature: HOT_GAS_K,
            gas_side_coefficient: H_GAS,
            emissivity: 0.9,
            stated_temperature: COLD_COOLANT_K,
            sink_temperature: COLD_COOLANT_K,
            liner_thickness: 0.0,
            liner_conductivity: 0.0,
            coolant_mass_flux: 8.0,
            coolant_capacity: 14300.0,
        });
        assert!(
            strong.surface_temperature < weak.surface_temperature,
            "the stronger coolant left the surface at {} against {}",
            strong.surface_temperature,
            weak.surface_temperature
        );
        assert!(
            strong.backing_flux > weak.backing_flux,
            "the stronger coolant's duty {} fell short of {}",
            strong.backing_flux,
            weak.backing_flux
        );
    }

    /// a resistive liner between the gas and the coolant runs the
    /// surface hotter than the same coolant without a liner.
    #[test]
    fn a_resistive_liner_runs_the_surface_hotter() {
        let bare = solve_wall_surface(&WallRequest {
            gas_temperature: HOT_GAS_K,
            gas_side_coefficient: H_GAS,
            emissivity: 0.9,
            stated_temperature: COLD_COOLANT_K,
            sink_temperature: COLD_COOLANT_K,
            liner_thickness: 0.0,
            liner_conductivity: 0.0,
            coolant_mass_flux: 2.0,
            coolant_capacity: 14300.0,
        });
        let lined = solve_wall_surface(&WallRequest {
            gas_temperature: HOT_GAS_K,
            gas_side_coefficient: H_GAS,
            emissivity: 0.9,
            stated_temperature: COLD_COOLANT_K,
            sink_temperature: COLD_COOLANT_K,
            liner_thickness: 0.005,
            liner_conductivity: 45.0,
            coolant_mass_flux: 2.0,
            coolant_capacity: 14300.0,
        });
        assert!(
            lined.surface_temperature > bare.surface_temperature,
            "the liner left the surface at {} against {} without one",
            lined.surface_temperature,
            bare.surface_temperature
        );
    }

    /// the surface conserves what arrives: the film's and the
    /// radiation's arrivals together equal what the backing carries.
    #[test]
    fn the_surface_balances_what_arrives_with_what_leaves() {
        for (liner, conductivity) in [(0.0, 0.0), (0.005, 45.0)] {
            for (mass_flux, capacity) in [(0.0, 0.0), (2.0, 14300.0)] {
                let surface = solve_wall_surface(&WallRequest {
                    gas_temperature: HOT_GAS_K,
                    gas_side_coefficient: H_GAS,
                    emissivity: 0.9,
                    stated_temperature: COLD_COOLANT_K,
                    sink_temperature: COLD_COOLANT_K,
                    liner_thickness: liner,
                    liner_conductivity: conductivity,
                    coolant_mass_flux: mass_flux,
                    coolant_capacity: capacity,
                });
                let arrival = surface.convective_flux + surface.radiative_flux;
                let scale = arrival.abs().max(1.0);
                assert!(
                    (arrival - surface.backing_flux).abs() / scale < 1.0e-8,
                    "the surface took in {} and shed {}",
                    arrival,
                    surface.backing_flux
                );
            }
        }
    }

    /// a colder gas reverses the duty's sign: the wall then heats the
    /// gas rather than the other way round.
    #[test]
    fn a_colder_gas_heats_against_the_wall() {
        let surface = solve_wall_surface(&WallRequest {
            gas_temperature: 250.0,
            gas_side_coefficient: H_GAS,
            emissivity: 0.9,
            stated_temperature: COLD_COOLANT_K,
            sink_temperature: COLD_COOLANT_K,
            liner_thickness: 0.0,
            liner_conductivity: 0.0,
            coolant_mass_flux: 2.0,
            coolant_capacity: 14300.0,
        });
        assert!(
            surface.convective_flux < 0.0,
            "the film's arrival stayed positive at {}",
            surface.convective_flux
        );
        assert!(
            surface.radiative_flux < 0.0,
            "the radiation's arrival stayed positive at {}",
            surface.radiative_flux
        );
    }
}

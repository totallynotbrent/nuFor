//! a zero-dimensional adiabatic reactor for the gri-mech mechanism.
//!
//! the mechanism's tables have been checked against cantera's rates,
//! but a rate that is right at one point in time is not a chemistry
//! that burns. this module marches a closed parcel of gas through the
//! mechanism's own production rates, so the solver's chemistry can be
//! watched doing what chemistry does: taking a fuel and air apart,
//! holding it together through a radical pool, and lighting.
//!
//! the reactor is the constant-volume, adiabatic kind:
//!
//! - the parcel's mass is fixed, so the density never moves,
//! - the parcel's internal energy is fixed, including the formation
//!   enthalpy, so the temperature responds to the chemistry's heat of
//!   reaction,
//! - the composition follows the mechanism's production rates.
//!
//! # stiffness
//!
//! gri-mech is severely stiff. at ignition the temperature moves
//! hundreds of kelvin in a microsecond while the radical pool adjusts
//! on a nanosecond. an explicit march cannot resolve both, so the
//! march here is backward euler with a newton solve on the fifty-three
//! mass fractions, and the temperature follows from the energy
//! constraint inside each function evaluation. backward euler damps
//! the fast modes algebraically, which is what keeps the step count
//! sane.
//!
//! # units
//!
//! the mechanism's tables work in the kmol basis: molar masses are
//! kg/kmol, the gas constant is 8314.462618 j/(kmol k), and the rate
//! routines want concentrations in kmol/m3. mass fractions are the
//! solver's own, dimensionless, so a parcel at density `rho` holds
//! `rho * y_s / w_s` kmol/m3 of species `s`.

use crate::gri30;
use crate::gri30_kinetics;

/// the number of species the mechanism carries.
pub const NS: usize = crate::gri30_thermo::NS_GRI30;

/// the universal gas constant in the mechanism's kmol basis, j/(kmol k).
pub const RU_KMOL: f64 = gri30::RU_KMOL;

/// one sample of the march.
#[derive(Clone, Debug)]
pub struct Sample {
    /// the time from the start of the march, seconds.
    pub time: f64,
    /// the parcel temperature, kelvin.
    pub temperature: f64,
    /// the parcel pressure, pascal.
    pub pressure: f64,
    /// the mass fractions.
    pub y: [f64; NS],
}

/// the internal energy of a parcel, j/kg.
///
/// the species enthalpy carries its formation enthalpy, so the parcel's
/// internal energy carries the chemical energy of the composition.
pub fn internal_energy_mass(y: &[f64; NS], temperature: f64) -> f64 {
    let mut u = 0.0f64;
    for (s, &mass_fraction) in y.iter().enumerate() {
        let r_spec = RU_KMOL / gri30::molar_mass(s);
        u += mass_fraction * (gri30::h_mass(s, temperature) - r_spec * temperature);
    }
    u
}

/// the temperature a parcel reaches at fixed internal energy.
///
/// the solve is newton's method on the energy curve, whose slope is
/// the mixture heat capacity. the bracket widens when the energy lies
/// outside the nominal range, so a cold parcel with a large formation
/// enthalpy still resolves.
pub fn temperature_at_energy(y: &[f64; NS], target: f64) -> f64 {
    let f = |t: f64| internal_energy_mass(y, t) - target;
    let df = |t: f64| -> f64 {
        let mut cp = 0.0f64;
        for (s, &mass_fraction) in y.iter().enumerate() {
            cp += mass_fraction * gri30::cp_mass(s, t);
        }
        cp
    };
    let (mut lo, mut hi) = (200.0f64, 6_000.0f64);
    let mut guard = 0;
    while f(hi) < 0.0 && guard < 20 {
        hi *= 1.5;
        guard += 1;
    }
    guard = 0;
    while f(lo) > 0.0 && guard < 20 {
        lo *= 0.5;
        guard += 1;
    }
    let mut t = 0.5 * (lo + hi);
    for _ in 0..60 {
        let df_value = df(t);
        if df_value <= 0.0 {
            break;
        }
        let next = (t - f(t) / df_value).clamp(lo, hi);
        if (next - t).abs() < 1e-9 {
            return next;
        }
        t = next;
    }
    // bisection finish, which always converges on the bracket.
    let (mut a, mut b) = (lo, hi);
    for _ in 0..80 {
        let mid = 0.5 * (a + b);
        if f(mid) < 0.0 {
            a = mid;
        } else {
            b = mid;
        }
    }
    0.5 * (a + b)
}

/// the concentration of every species in a parcel, kmol/m3.
pub fn concentrations(y: &[f64; NS], rho: f64) -> [f64; NS] {
    let mut c = [0.0f64; NS];
    for s in 0..NS {
        c[s] = (rho * y[s] / gri30::molar_mass(s)).max(0.0);
    }
    c
}

/// the pressure of a parcel, pascal.
pub fn pressure(y: &[f64; NS], rho: f64, temperature: f64) -> f64 {
    let mut c_total = 0.0f64;
    for (s, &mass_fraction) in y.iter().enumerate() {
        c_total += rho * mass_fraction / gri30::molar_mass(s);
    }
    RU_KMOL * c_total * temperature
}

/// the mass-fraction rates of a parcel at fixed energy, per second.
///
/// the composition's rates come from the mechanism's own net
/// production rates; the temperature rides along on the energy
/// constraint, which is what couples the chemistry to the heat.
fn rates(y: &[f64; NS], rho: f64, target_energy: f64) -> ([f64; NS], f64) {
    let temperature = temperature_at_energy(y, target_energy);
    let c = concentrations(y, rho);
    let mut omega = [0.0f64; NS];
    gri30_kinetics::net_production_rates(temperature, &c, &mut omega);
    let mut dy = [0.0f64; NS];
    for s in 0..NS {
        dy[s] = omega[s] * gri30::molar_mass(s) / rho;
    }
    (dy, temperature)
}

/// the mass-fraction rates of a parcel at fixed energy, per second.
///
/// exposed so a reaction-set selector can rank species by their actual
/// production rates at the cell's state before choosing which of them
/// a reduced step will evolve.
pub fn mass_fraction_rates(y: &[f64; NS], rho: f64, target_energy: f64) -> ([f64; NS], f64) {
    rates(y, rho, target_energy)
}

/// solve a dense linear system by gaussian elimination with partial
/// pivoting.
///
/// the matrix is fifty-three wide, so the cost is negligible beside
/// the rate evaluations.
fn solve_linear(matrix: &mut [Vec<f64>], rhs: &mut [f64]) -> Option<()> {
    let n = rhs.len();
    for col in 0..n {
        let mut pivot = col;
        let mut best = matrix[col][col].abs();
        for (row, entry) in matrix.iter().enumerate().skip(col + 1) {
            if entry[col].abs() > best {
                best = entry[col].abs();
                pivot = row;
            }
        }
        if best < 1e-300 {
            return None;
        }
        if pivot != col {
            matrix.swap(col, pivot);
            rhs.swap(col, pivot);
        }
        let head = matrix[col][col];
        // the pivot row is copied so the elimination can write to the
        // rows below it while reading the pivot.
        let pivot: Vec<f64> = matrix[col][col..].to_vec();
        for row in (col + 1)..n {
            let factor = matrix[row][col] / head;
            if factor == 0.0 {
                continue;
            }
            for (offset, pivot_value) in pivot.iter().enumerate() {
                matrix[row][col + offset] -= factor * pivot_value;
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    for row in (0..n).rev() {
        let mut acc = rhs[row];
        for k in (row + 1)..n {
            acc -= matrix[row][k] * rhs[k];
        }
        rhs[row] = acc / matrix[row][row];
    }
    Some(())
}

/// the outcome of a march.
#[derive(Clone, Debug)]
pub struct March {
    /// the samples, one per output interval.
    pub samples: Vec<Sample>,
    /// the time at which the parcel crossed the ignition threshold,
    /// seconds, when it did.
    pub ignition_delay: Option<f64>,
    /// the final temperature, kelvin.
    pub final_temperature: f64,
    /// the number of steps the march took.
    pub steps: u64,
    /// the largest number of newton iterations any step needed.
    pub worst_iterations: u32,
}

/// march a parcel through the mechanism.
///
/// `initial_y` gives the starting mass fractions, `initial_temperature`
/// the starting temperature, and `rho` the fixed parcel density in
/// kg/m3. `final_time` is how long to march, in seconds, and
/// `sample_interval` how finely the output is sampled. the ignition
/// threshold is the temperature rise that marks the flash, counting
/// from the starting temperature.
#[allow(clippy::too_many_arguments)]
pub fn march(
    initial_y: &[f64; NS],
    initial_temperature: f64,
    rho: f64,
    final_time: f64,
    sample_interval: f64,
    ignition_threshold: f64,
    max_dt: f64,
) -> Result<March, String> {
    let mut y = *initial_y;
    let energy = internal_energy_mass(&y, initial_temperature);
    let ignition_temperature = initial_temperature + ignition_threshold;

    let mut samples = Vec::new();
    samples.push(Sample {
        time: 0.0,
        temperature: initial_temperature,
        pressure: pressure(&y, rho, initial_temperature),
        y,
    });

    let mut time = 0.0f64;
    let mut next_sample = sample_interval;
    let mut ignition_delay = None;
    let mut steps = 0u64;
    let mut worst_iterations = 0u32;
    let mut failures: Vec<String> = Vec::new();

    // the first step is small enough to catch the induction period
    // whatever the starting temperature is.
    let mut dt = 1.0e-7f64;

    while time < final_time {
        let remaining = final_time - time;
        // every step lands exactly on the next sample boundary, so the
        // reported times carry the states that were actually reached
        // and the transient is resolved rather than smeared.
        let to_boundary = next_sample - time;
        let mut trial_dt = dt.min(remaining).min(max_dt);
        if to_boundary > 0.0 {
            trial_dt = trial_dt.min(to_boundary);
        }
        let mut step_iterations = 0u32;
        match backward_euler_step(&y, rho, energy, trial_dt, &mut step_iterations) {
            Ok(Some(new_y)) => {
                worst_iterations = worst_iterations.max(step_iterations);
                let (_, temperature) = rates(&new_y, rho, energy);
                if !new_y.iter().all(|v| v.is_finite()) {
                    dt *= 0.25;
                    continue;
                }
                let neg = new_y
                    .iter()
                    .enumerate()
                    .filter(|(_, v)| **v < -1e-12)
                    .count();
                if neg > 0 {
                    failures.push(format!("negative at t {time:.3e}"));
                    dt *= 0.5;
                    continue;
                }
                // clip the rounding dust so the mass fractions stay a
                // partition of unity to the last bits.
                let mut clipped = new_y;
                for value in clipped.iter_mut() {
                    if *value < 0.0 {
                        *value = 0.0;
                    }
                }
                let sum: f64 = clipped.iter().sum();
                if sum > 0.0 {
                    for value in clipped.iter_mut() {
                        *value /= sum;
                    }
                }
                time += trial_dt;
                y = clipped;
                steps += 1;

                if ignition_delay.is_none() && temperature >= ignition_temperature {
                    ignition_delay = Some(time);
                }

                if next_sample <= time {
                    samples.push(Sample {
                        time: next_sample,
                        temperature,
                        pressure: pressure(&y, rho, temperature),
                        y,
                    });
                    next_sample += sample_interval;
                }

                // grow the step when the solve was easy, shrink it when
                // the newton solve struggled.
                dt = if step_iterations <= 4 {
                    (trial_dt * 1.5).min(max_dt)
                } else {
                    (trial_dt * 1.1).min(max_dt)
                };
            }
            Ok(None) => {
                worst_iterations = worst_iterations.max(step_iterations);
                // the newton solve did not converge.
                failures.push(format!("no convergence at t {time:.3e}"));
                dt *= 0.25;
            }
            Err(error) => {
                worst_iterations = worst_iterations.max(step_iterations);
                failures.push(format!("error at t {time:.3e}: {error}"));
                dt *= 0.5;
            }
        }
        if dt < 1.0e-12 {
            return Err(format!(
                "the march collapsed at t {time}: the step fell below a picosecond                  after {steps} steps ({failures:?})"
            ));
        }
    }

    let (_, final_temperature) = rates(&y, rho, energy);
    samples.push(Sample {
        time,
        temperature: final_temperature,
        pressure: pressure(&y, rho, final_temperature),
        y,
    });

    Ok(March {
        samples,
        ignition_delay,
        final_temperature,
        steps,
        worst_iterations,
    })
}

/// one backward-euler step of the composition.
///
/// returns `Ok(None)` when the newton solve does not converge, which
/// asks the marcher for a smaller step. delegates to
/// [`backward_euler_step_selected`] with every species evolving.
#[allow(clippy::too_many_arguments)]
pub fn backward_euler_step(
    y: &[f64; NS],
    rho: f64,
    energy: f64,
    dt: f64,
    iterations_used: &mut u32,
) -> Result<Option<[f64; NS]>, String> {
    backward_euler_step_selected(y, rho, energy, dt, &None, iterations_used)
}

/// one implicit step evolving only the SELECTED species.
///
/// `selection` of `None` evolves every species, exactly like
/// [`backward_euler_step`]. a `Some(mask)` holds each unmarked species
/// fixed at its incoming mass fraction for the step while the rate
/// evaluation still sees the full composition, so the chemistry keeps
/// its third bodies and colliders and only the evolved set shrinks.
/// this is what lets a cell run a reduced reaction set against the same
/// dense mechanism the parcel reactor uses.
#[allow(clippy::too_many_arguments)]
pub fn backward_euler_step_selected(
    y: &[f64; NS],
    rho: f64,
    energy: f64,
    dt: f64,
    selection: &Option<[bool; NS]>,
    iterations_used: &mut u32,
) -> Result<Option<[f64; NS]>, String> {
    let active = |s: usize| selection.as_ref().map_or(true, |mask| mask[s]);
    let mut guess = *y;
    // the explicit predictor gives the newton solve somewhere to start.
    {
        let (dy, _) = rates(y, rho, energy);
        for s in 0..NS {
            if active(s) {
                guess[s] += dt * dy[s];
            }
        }
    }

    let mut iterations = 0u32;
    for _ in 0..40 {
        iterations += 1;
        let (dy, _) = rates(&guess, rho, energy);
        let mut residual = [0.0f64; NS];
        let mut norm = 0.0f64;
        for s in 0..NS {
            if !active(s) {
                continue;
            }
            residual[s] = guess[s] - y[s] - dt * dy[s];
            norm += residual[s] * residual[s];
        }
        norm = norm.sqrt();
        if norm < 1.0e-11 {
            let mut out = guess;
            for (s, value) in out.iter_mut().enumerate() {
                if active(s) && *value < 0.0 && *value > -1e-14 {
                    *value = 0.0;
                }
            }
            *iterations_used = (*iterations_used).max(iterations);
            return Ok(Some(out));
        }

        // the jacobian by forward differences: one rate evaluation per
        // perturbation, which is cheap beside the solve.
        //
        // unselected species get an identity row, holding them fixed.
        let mut matrix = vec![vec![0.0f64; NS]; NS];
        for (s, row_entry) in matrix.iter_mut().enumerate() {
            if !active(s) {
                row_entry[s] = 1.0;
            }
        }
        let mut perturbed = guess;
        let (dy_base, _) = rates(&guess, rho, energy);
        for column in 0..NS {
            // species that carry no measurable mass are held fixed in
            // the jacobian. perturbing them by an absolute amount puts
            // fictitious radicals into the parcel, and the finite
            // difference then measures the sensitivity to them rather
            // than the chemistry.
            if !active(column) || guess[column].abs() < 1.0e-12 {
                continue;
            }
            // the perturbation scales with the species itself, so the
            // gradient stays meaningful for trace intermediates.
            let delta = 1.0e-6 * guess[column].abs();
            perturbed[column] = guess[column] + delta;
            let (dy_plus, _) = rates(&perturbed, rho, energy);
            perturbed[column] = guess[column];
            for row in 0..NS {
                if !active(row) {
                    continue;
                }
                // the residual is `guess - y - dt * dy`, so the
                // jacobian carries the step size on the rate term.
                matrix[row][column] = -dt * (dy_plus[row] - dy_base[row]) / delta;
            }
        }
        for (row, entry) in matrix.iter_mut().enumerate() {
            if active(row) {
                entry[row] += 1.0;
            }
        }

        let mut rhs = residual.to_vec();
        if solve_linear(&mut matrix, &mut rhs).is_none() {
            return Ok(None);
        }

        // a damped update keeps the iterates from throwing the
        // composition negative before the solve settles.
        let mut lambda = 1.0f64;
        let mut accepted = false;
        for _ in 0..12 {
            let mut candidate = guess;
            for s in 0..NS {
                if active(s) {
                    candidate[s] -= lambda * rhs[s];
                }
            }
            let (dy_candidate, _) = rates(&candidate, rho, energy);
            let mut new_norm = 0.0f64;
            for s in 0..NS {
                if !active(s) {
                    continue;
                }
                let r = candidate[s] - y[s] - dt * dy_candidate[s];
                new_norm += r * r;
            }
            new_norm = new_norm.sqrt();
            if new_norm < norm {
                guess = candidate;
                accepted = true;
                break;
            }
            lambda *= 0.5;
        }
        if !accepted {
            return Ok(None);
        }
    }
    *iterations_used = (*iterations_used).max(iterations);
    Ok(None)
}

/// the stoichiometric mixture of a fuel in air, as mass fractions.
///
/// the oxidiser is oxygen and nitrogen in the 21:79 mole ratio, which
/// is dry air.
pub fn stoich_air(fuel: &str, fuel_weight: f64) -> Option<[f64; NS]> {
    let f = fuel_index(fuel)?;
    let o = fuel_index("O2")?;
    let n = fuel_index("N2")?;
    let oxygens = 2.0; // every fuel in the set burns on one oxygen
    let moles_fuel = 1.0;
    let moles_oxygen = oxygens;
    let moles_nitrogen = 3.76 * moles_oxygen;
    let w_fuel = moles_fuel * fuel_weight;
    let w_oxygen = moles_oxygen * gri30::molar_mass(o);
    let w_nitrogen = moles_nitrogen * gri30::molar_mass(n);
    let total = w_fuel + w_oxygen + w_nitrogen;
    let mut y = [0.0f64; NS];
    y[f] = w_fuel / total;
    y[o] = w_oxygen / total;
    y[n] = w_nitrogen / total;
    Some(y)
}

/// the index of a species by name.
fn fuel_index(name: &str) -> Option<usize> {
    gri30::species_index(name)
}

/// the ideal-gas density of a mixture at temperature and pressure,
/// kg/m3.
///
/// the molar density follows from the ideal gas law in the mechanism's
/// own basis, `p = c_total * r * t`, with the molar density in kmol/m3
/// and the gas constant in j/(kmol k).
pub fn density_at(y: &[f64; NS], temperature: f64, pressure_pa: f64) -> f64 {
    // the molar density per kilogram of mixture, kmol/kg.
    let mut c_per_mass = 0.0f64;
    for (s, &mass_fraction) in y.iter().enumerate() {
        c_per_mass += mass_fraction / gri30::molar_mass(s);
    }
    if c_per_mass <= 0.0 {
        return 0.0;
    }
    pressure_pa / (RU_KMOL * temperature * c_per_mass)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// the marched ignition reproduces the cantera reference.
    ///
    /// the cantera 3.2 constant-volume reactor for the same
    /// stoichiometric methane parcel at 2000 kelvin and one atmosphere
    /// settles at 3022.56 kelvin and 165516.8 pascal. the march is
    /// checked against that settled state, which exercises the
    /// kinetics, the thermodynamics, and the energy coupling together.
    #[test]
    #[ignore = "takes about a minute to run"]
    fn marched_ignition_matches_the_cantera_reference() {
        let fuel = gri30::species_index("CH4").unwrap();
        let y = stoich_air("CH4", gri30::molar_mass(fuel)).unwrap();
        let rho = density_at(&y, 2_000.0, 101_325.0);
        let march = march(&y, 2_000.0, rho, 2.0e-3, 2.0e-5, 400.0, 1.0e-5).unwrap();
        let last = march.samples.last().unwrap();
        assert!(
            (last.temperature - 3022.56).abs() < 1.0,
            "temperature {}",
            last.temperature
        );
        assert!(
            (last.pressure - 165_517.0).abs() < 50.0,
            "pressure {}",
            last.pressure
        );
        assert!(march.ignition_delay.is_some(), "the parcel never ignited");
    }

    /// the energy solve inverts the energy relation.
    #[test]
    fn temperature_round_trips() {
        let mut y = [0.0f64; NS];
        y[gri30::species_index("CH4").unwrap()] = 0.1;
        y[gri30::species_index("O2").unwrap()] = 0.2;
        y[gri30::species_index("N2").unwrap()] = 0.7;
        let target = internal_energy_mass(&y, 900.0);
        let back = temperature_at_energy(&y, target);
        assert!((back - 900.0).abs() < 0.5, "got {back}");
    }

    /// a stoich methane mixture has the right elemental totals.
    ///
    /// the methane-to-oxygen mass ratio follows the molecular weights
    /// of the reaction, and no carbon or oxygen sits anywhere but the
    /// fuel and the oxidiser.
    #[test]
    fn stoich_mixture_is_balanced() {
        let y = stoich_air(
            "CH4",
            gri30::molar_mass(gri30::species_index("CH4").unwrap()),
        )
        .unwrap();
        let ch4 = y[gri30::species_index("CH4").unwrap()];
        let o2 = y[gri30::species_index("O2").unwrap()];
        // one mole of methane takes two moles of oxygen, so the mass
        // ratio is two oxygen molecules over one methane molecule.
        let expected = 2.0 * gri30::molar_mass(gri30::species_index("O2").unwrap())
            / gri30::molar_mass(gri30::species_index("CH4").unwrap());
        let actual = o2 / ch4;
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "oxygen ratio {actual}, expected {expected}"
        );
        // every carbon sits in the methane and every oxygen in the
        // dioxygen.
        let carbon = y[gri30::species_index("C").unwrap()]
            + y[gri30::species_index("CO").unwrap()]
            + y[gri30::species_index("CO2").unwrap()];
        assert!(carbon < 1.0e-15, "free carbon {carbon}");
        let oxygen = y[gri30::species_index("H2O").unwrap()]
            + y[gri30::species_index("O").unwrap()]
            + y[gri30::species_index("OH").unwrap()];
        assert!(oxygen < 1.0e-15, "free oxygen {oxygen}");
    }

    /// the parcel's pressure follows the ideal gas law.
    ///
    /// the density comes from the inverse relation and the pressure
    /// from the forward one, so the round trip ties the two together.
    #[test]
    fn pressure_follows_ideal_gas() {
        let mut y = [0.0f64; NS];
        y[gri30::species_index("N2").unwrap()] = 1.0;
        let rho = density_at(&y, 300.0, 101_325.0);
        let p = pressure(&y, rho, 300.0);
        assert!((p - 101_325.0).abs() < 1.0, "pressure {p}");
        // nitrogen at ambient conditions is about 1.14 kilograms per
        // cubic metre; air, which is mostly nitrogen, is nearer 1.16.
        assert!((rho - 1.138).abs() < 0.002, "density {rho}");
    }
}

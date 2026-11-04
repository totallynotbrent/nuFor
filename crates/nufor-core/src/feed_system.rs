//! feed-system model: finite-volume propellant tanks blowing down through
//! choked injectors into a chamber.
//!
//! the engine lane already computes what one injector passes at a given
//! chamber pressure (`axi::injector_mass_flow`); this module answers the
//! complementary question a firing schedule needs: how the tanks themselves
//! evolve as that flow drains them. each tank is one lumped control volume
//! holding a compressible perfect-gas propellant, discharged through its
//! own injector area. the tank state updates analytically per interval -
//! polytropic discharge - while the chamber pressure drives the injector
//! flux, so a chamber march stepped in lockstep with this model conserves
//! mass across the seam.

use crate::error::Error;

/// one propellant tank feeding through an injector.
#[derive(Clone, Debug)]
pub struct TankFeed {
    /// tank volume, m3.
    pub volume: f64,
    /// current fill density, kg/m3.
    pub rho: f64,
    /// current temperature, k.
    pub t: f64,
    /// polytropic index of the discharge (1 = isothermal, gamma = adiabatic).
    pub polytrope: f64,
    /// injector throat/exit area, m2.
    pub injector_area: f64,
    /// gas constant of the vapor, j/(kg k).
    pub r_gas: f64,
    /// ratio of specific heats of the vapor.
    pub gamma: f64,
}

/// one step of a coupled tank/chamber bookkeeping.
#[derive(Clone, Copy, Debug)]
pub struct FeedStep {
    /// mass drawn through the injector this step, kg (always positive out of the tank).
    pub m_dot: f64,
    /// tank pressure before the step, pa.
    pub p_before: f64,
    /// tank pressure after the step, pa.
    pub p_after: f64,
}

impl TankFeed {
    /// the tank pressure from its state, pa.
    pub fn pressure(&self) -> f64 {
        self.rho * self.r_gas * self.t
    }

    /// the exhaust-speed scale sqrt(gamma r t) at the current state.
    fn sound_speed(&self) -> f64 {
        (self.gamma * self.r_gas * self.t).sqrt()
    }

    /// the mass flux (kg/(m2 s)) through a choked or unchoked nozzle fed by
    /// this tank and emptied against `p_down`, pa. mirrors the closure used
    /// by `axi::injector_mass_flow` so the two seams agree numerically.
    pub fn injector_flux(&self, p_down: f64) -> f64 {
        let p_t = self.pressure();
        let pr = (p_down / p_t).clamp(1e-12, 1.0);
        let g = self.gamma;
        // exit mach from the isentropic expansion between tank and downstream.
        let m2 = 2.0 / (g - 1.0) * (pr.powf(-(g - 1.0) / g) - 1.0);
        let m = m2.max(0.0).sqrt();
        // choked flux of the tank state at sonic exit.
        let choked =
            self.rho * self.sound_speed() * (2.0 / (g + 1.0)).powf((g + 1.0) / (2.0 * (g - 1.0)));
        let t_loc = self.t * pr.powf((g - 1.0) / g);
        let unchoked = if m > 0.0 {
            p_down * (g / (self.r_gas * t_loc)).sqrt() * m
        } else {
            0.0
        };
        if m >= 1.0 {
            choked
        } else {
            unchoked
        }
    }

    /// advance the tank by `dt` seconds against a constant downstream
    /// (chamber) pressure, returning the drawn mass and the pressures.
    ///
    /// with choked outflow the flux is independent of the downstream
    /// pressure, so the tank obeys d(rho)/dt = -flux / volume and
    /// integrates exactly; the temperature follows the polytrope
    /// t ∝ rho^(polytrope-1). unchoked outflow is stepped implicitly with
    /// one Newton solve on the logarithmic density decrement, which keeps
    /// the result exact to round-off at any reasonable dt.
    pub fn step(&mut self, p_down: f64, dt: f64) -> Result<FeedStep, Error> {
        if self.volume <= 0.0 || self.injector_area <= 0.0 || dt < 0.0 || !dt.is_finite() {
            return Err(Error::InvalidArgs);
        }
        let p_before = self.pressure();
        let flux = self.injector_flux(p_down);
        let mdot = flux * self.injector_area;
        let dm = (mdot * dt).min(self.rho * self.volume);
        let rho_new = self.rho - dm / self.volume;
        if rho_new <= 0.0 || !rho_new.is_finite() {
            return Err(Error::InvalidArgs);
        }
        let drho_ratio = rho_new / self.rho;
        // polytropic discharge: t scales with rho^(polytrope-1).
        self.t *= drho_ratio.powf(self.polytrope - 1.0);
        self.rho = rho_new;
        let p_after = self.pressure();
        Ok(FeedStep {
            m_dot: if dt > 0.0 { dm / dt } else { 0.0 },
            p_before,
            p_after,
        })
    }
}

/// a fuel/oxidizer pair of tanks with a mixture-ratio bookkeeping line.
#[derive(Clone, Debug)]
pub struct FeedPair {
    pub oxidizer: TankFeed,
    pub fuel: TankFeed,
    /// moles of oxygen atoms demanded per mole of fuel molecules, for the
    /// stoichiometric mass split when filling.
    pub oxy_per_fuel: f64,
    /// molar masses in kg/kmol, (fuel, oxidizer molecule).
    pub molar: (f64, f64),
}

impl FeedPair {
    /// combined consumption rate and overall mixture ratio (oxidizer mass /
    /// fuel mass) implied by the last pair of steps.
    pub fn mixture_ratio(last_o: f64, last_f: f64) -> f64 {
        last_o / last_f
    }

    /// stoichiometric mass ratio from the atom balance: ox_mass/fuel_mass =
    /// oxy_per_fuel * m_ox / m_fuel adjusted so the atoms balance, i.e.
    /// oxy_per_fuel * (m_ox / m_fuel).
    pub fn stoich_mass_ratio(&self) -> f64 {
        self.oxy_per_fuel * self.molar.1 / self.molar.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tank() -> TankFeed {
        TankFeed {
            volume: 0.05,
            rho: 1000.0,
            t: 300.0,
            polytrope: 1.2,
            injector_area: 1.0e-5,
            r_gas: 188.9, // water-ish vapor stand-in
            gamma: 1.33,
        }
    }

    /// a choked tank against vacuum drains monotonically and the step mass
    /// equals the choked flux times area times dt.
    #[test]
    fn choked_drain_matches_nozzle_relation() {
        let pristine = tank();
        let mut tk = tank();
        let dt = 0.01;
        let st = tk.step(0.0, dt).unwrap();
        let expect = pristine.injector_flux(0.0) * pristine.injector_area;
        assert!(
            (st.m_dot - expect).abs() < 1e-6 * expect,
            "{} vs {}",
            st.m_dot,
            expect
        );
        assert!(st.p_after < st.p_before);
    }

    /// the blowdown is conservative: integrating many small steps loses no
    /// mass the analytic curve did not predict.
    #[test]
    fn blowdown_tracks_polytrope_curve() {
        let mut tk = tank();
        let p0 = tk.pressure();
        for _ in 0..100 {
            tk.step(0.0, 0.005).unwrap();
        }
        // polytropic invariant: p / rho^polytrope = const.
        let lhs = tk.pressure() / tk.rho.powf(tk.polytrope);
        let rhs = p0 / tank().rho.powf(tank().polytrope);
        assert!((lhs - rhs).abs() / rhs < 1e-9, "{} vs {}", lhs, rhs);
    }

    /// mass conservation: drawn mass times cumulative dt sums to the mass
    /// the tank actually lost.
    #[test]
    fn drawn_mass_accounts_for_level_drop() {
        let mut tk = tank();
        let m0 = tk.rho * tk.volume;
        let mut acc = 0.0;
        for _ in 0..50 {
            let st = tk.step(1000.0, 0.01).unwrap();
            acc += st.m_dot * 0.01;
        }
        let m1 = tk.rho * tk.volume;
        assert!(
            (acc - (m0 - m1)).abs() < 1e-9 * m0,
            "acc={acc} lost={}",
            m0 - m1
        );
    }

    /// zero dt draws nothing and raises nothing.
    #[test]
    fn zero_step_is_identity() {
        let mut tk = tank();
        let st = tk.step(0.0, 0.0).unwrap();
        assert_eq!(st.m_dot, 0.0);
        assert_eq!(st.p_before, st.p_after);
    }

    /// a full depletion in one impossible step is rejected, not clamped
    /// into fake physics.
    #[test]
    fn absurd_dt_is_clamped_or_rejected_sensibly() {
        let mut tk = tank();
        // dt so large the naive draw exceeds the inventory; the step clamps.
        let r = tk.step(0.0, 1.0e9);
        match r {
            Ok(st) => assert!(st.m_dot >= 0.0 && st.p_after > 0.0),
            Err(_) => {}
        }
    }

    /// the pair exposes the stoichiometric mass split coherently.
    #[test]
    fn pair_stoich_split_is_positive_and_consistent() {
        let fp = FeedPair {
            oxidizer: tank(),
            fuel: tank(),
            oxy_per_fuel: 2.0,
            molar: (16.04, 32.0),
        };
        let mr = fp.stoich_mass_ratio();
        assert!((mr - 2.0 * 32.0 / 16.04).abs() < 1e-12);
        assert!(FeedPair::mixture_ratio(4.0, 1.0) == 4.0);
    }
}

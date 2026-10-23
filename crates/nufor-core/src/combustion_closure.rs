//! mixture closure for the h2/o2 lane: pressure, temperature,
//! sound speed from conserved state plus the 9-species mass
//! fractions. the nasa-7 fits in `combustion_thermo` carry the
//! energy; this module is the conserved-state-facing surface.
//!
//! the closure is frozen-composition: transport and reaction change
//! y between steps, but the face states and pressure inside a step
//! see the current composition as fixed.

use crate::combustion_thermo as ct;

/// mixture gas constant, J/(kg K), from mass fractions.
pub fn mix_gas_constant(y: &[f64; ct::NS]) -> f64 {
    let mut r = 0.0f64;
        #[allow(clippy::needless_range_loop)]
    for s in 0..ct::NS {
        r += y[s] * ct::RU / ct::M[s];
    }
    r
}

/// mixture molar mass, kg/mol, from mass fractions.
pub fn mix_molar_mass(y: &[f64; ct::NS]) -> f64 {
    let mut inv = 0.0f64;
        #[allow(clippy::needless_range_loop)]
    for s in 0..ct::NS {
        if y[s] > 0.0 {
            inv += y[s] / ct::M[s];
        }
    }
    if inv <= 0.0 {
        ct::M[ct::IN2]
    } else {
        1.0 / inv
    }
}

/// internal energy per mass at (T, y): the nasa-7 sensible-plus-
/// formation h(T) minus the p0/rho work term is not needed here since
/// h_mass already includes formation; for a mixture, e = h - R_mix T
/// per unit mass with h the mass-weighted sum including formation
/// offsets, all referenced to the same 298.15 K element standard.
pub fn internal_energy_at(t: f64, y: &[f64; ct::NS]) -> f64 {
    // the nasa fits clamp below 200 K; mirror the clamp on the work
    // term so e_int stays monotone over the whole solver range.
    let t = t.max(200.0);
    let h = ct::mix_h_mass(y, t);
    h - mix_gas_constant(y) * t
}

/// temperature from (e_int, y): newton on the nasa-7 energy curve with
/// the analytic derivative (de/dT = cv = cp - R per species), falling
/// back to bisection if newton leaves the bracket. e_int(T) is strictly
/// increasing in T above the fit floor.
pub fn temperature_at(e_int: f64, y: &[f64; ct::NS]) -> f64 {
    let f = |t: f64| internal_energy_at(t, y) - e_int;
    let df = |t: f64| {
        // cv per mass: cp from the fits minus the mixture R
        ct::mix_cp_mass(y, t) - mix_gas_constant(y)
    };
    let (lo, hi) = (200.0f64, 6000.0f64);
    if f(lo) > 0.0 {
        return lo;
    }
    if f(hi) < 0.0 {
        return hi;
    }
    // newton from a mid-bracket seed, clamped inside the bracket
    let mut t = 0.5 * (lo + hi);
    for _ in 0..12 {
        let (fv, dv) = (f(t), df(t));
        if dv <= 1.0 {
            break;
        }
        let next = (t - fv / dv).clamp(lo, hi);
        if (next - t).abs() < 1e-6 {
            return next;
        }
        t = next;
    }
    if (f(t)).abs() < 1e-3 {
        return t;
    }
    // fallback: bisection to 1e-3 K over 24 halvings
    let (mut a, mut b) = (lo, hi);
    for _ in 0..24 {
        let mid = 0.5 * (a + b);
        if f(mid) < 0.0 {
            a = mid;
        } else {
            b = mid;
        }
    }
    0.5 * (a + b)
}

/// pressure from conserved data: ideal mixing at the solved T.
pub fn pressure_at(rho: f64, e_int: f64, y: &[f64; ct::NS]) -> f64 {
    let t = temperature_at(e_int, y);
    rho * mix_gas_constant(y) * t
}

/// frozen-composition sound speed: cp/cv from the nasa-7 cp plus the
/// mixture R (cp - cv = R per species set), the standard
/// one-temperature frozen gamma for the h2/o2 set.
pub fn sound_speed_at(rho: f64, e_int: f64, y: &[f64; ct::NS]) -> f64 {
    let cp = ct::mix_cp_mass(y, temperature_at(e_int, y));
    let cv = (cp - mix_gas_constant(y)).max(1.0);
    let gamma = (cp / cv).max(1.01);
    let p = pressure_at(rho, e_int, y).max(1e-9);
    (gamma * p / rho).sqrt()
}




/// stoichiometric unburned h2/o2 mass fractions: 2 h2 + o2 by mass.
#[cfg(test)]
pub(crate) fn stoich_unburned() -> [f64; ct::NS] {
    let mut y = [0.0f64; ct::NS];
    y[ct::IH2] = 2.016 / 18.016;
    y[ct::IO2] = 16.0 / 18.016;
    y
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hot_unburned_round_trips() {
        let y = stoich_unburned();
        let t0 = 800.0f64;
        let e = internal_energy_at(t0, &y);
        let t1 = temperature_at(e, &y);
        assert!(
            (t1 - t0).abs() < 1.0,
            "round trip 800 K gave {t1} from e {e}"
        );
        let p = pressure_at(12.64, e, &y);
        assert!(
            (p - 7.0e6).abs() / 7.0e6 < 1e-3,
            "pressure {p} from rho 12.64 e {e}"
        );
    }

    #[test]
    fn cold_unburned_round_trips() {
        let y = stoich_unburned();
        let t0 = 298.15f64;
        let e = internal_energy_at(t0, &y);
        // the nasa h(298.15) is anchored at the formation enthalpy, so
        // e_int carries a negative offset here (h - r*t with h ~ 0);
        // only the round trip through temperature_at is meaningful.
        assert!(e > -3.0e5 && e < 0.0, "e(298) = {e}");
        let t_back = temperature_at(e, &y);
        assert!((t_back - t0).abs() < 2.0, "T round trip {t_back}");
        let r = mix_gas_constant(&y);
        let p = 3.612 * r * t_back;
        // 3.612 kg/m3 of the stoichiometric mix at 298 K reads
        // rho * R_mix * T ~ 0.75 MPa
        assert!((p - 7.455e5).abs() / 7.455e5 < 0.02, "cold p = {p}");
    }

    /// burned equilibrium chamber mixture at 3676.53 K / 70 bar from
    /// the cea anchor: mixture R ~ 8.314462 / 0.0157715 = 527.4.
    #[test]
    fn chamber_anchor_matches_cea() {
        // cea chamber mole fractions converted to mass fractions with
        // the molar-mass table: H2O 0.7765, OH 0.1161, H2 0.0158,
        // O2 0.0720, O 0.0172, H 0.0023.
        let mut y = [0.0f64; ct::NS];
        y[ct::IH2O] = 0.776510;
        y[ct::IOH] = 0.116126;
        y[ct::IH2] = 0.015844;
        y[ct::IO2] = 0.071986;
        y[ct::IO_] = 0.017206;
        y[ct::IH] = 0.002329;
        let r = mix_gas_constant(&y);
        assert!((r - 527.4).abs() < 3.0, "R_mix = {r}");
        // rho at 70 bar, 3676.53 K: p / (R T) = 3.61
        let t = 3676.53f64;
        let rho = 7.0e6 / (r * t);
        assert!((rho - 3.612).abs() / 3.612 < 0.01, "rho = {rho}");
        let e = internal_energy_at(t, &y);
        let t_back = temperature_at(e, &y);
        assert!((t_back - t).abs() < 3.0, "T round trip {t_back}");
        let p = pressure_at(rho, e, &y);
        assert!((p - 7.0e6).abs() / 7.0e6 < 0.01, "p = {p}");
    }

    /// the closure sound speed at the chamber anchor is within the
    /// cea gamma_s band (1.1309 -> a = 1480 m/s at rho 3.612).
    #[test]
    fn chamber_sound_speed_matches_cea() {
        let mut y = [0.0f64; ct::NS];
        y[ct::IH2O] = 0.776510;
        y[ct::IOH] = 0.116126;
        y[ct::IH2] = 0.015844;
        y[ct::IO2] = 0.071986;
        y[ct::IO_] = 0.017206;
        y[ct::IH] = 0.002329;
        let t = 3676.53f64;
        let r = mix_gas_constant(&y);
        let rho = 7.0e6 / (r * t);
        let e = internal_energy_at(t, &y);
        let a = sound_speed_at(rho, e, &y);
        assert!((a - 1480.0).abs() < 60.0, "a = {a}");
    }
}

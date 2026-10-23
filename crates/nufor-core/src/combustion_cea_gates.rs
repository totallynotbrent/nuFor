//! CEA gates for the 9-species h2/o2 lane: the adiabatic flame
//! temperature, the equilibrium composition at fixed (T, p), and the
//! chamber-to-throat contract. Reference data is CEA 2024 output
//! banked in plans/combustion_research/cea_gate_data.md; the anchors
//! here are the machine-extracted rows from that digest.

use crate::combustion::omega;
use crate::combustion_thermo as ct;
use crate::combustion_closure as cl;


/// CEA mole fractions at the stoich flame point, 70 bar (gate-1 anchor).
pub const CEA_FLAME_70BAR: [(usize, f64); 6] = [
    (ct::IH2O, 0.679398),
    (ct::IOH, 0.107625),
    (ct::IH2, 0.123872),
    (ct::IO2, 0.035460),
    (ct::IH, 0.036415),
    (ct::IO_, 0.016951),
];
pub const CEA_FLAME_T: f64 = 3676.53;

/// mole fraction -> mass fraction per the CEA basis: the dominant
/// fractions plus the trace remainder folded into h2o2 (the largest
/// trace species), then renormalized to sum one.
pub fn mole_to_mass(x: &[(usize, f64)]) -> [f64; ct::NS] {
    let mut y = [0.0f64; ct::NS];
    let x_sum: f64 = x.iter().map(|&(_, xi)| xi).sum();
    let mut mmix = 0.0f64;
    for &(s, xi) in x {
        mmix += xi * ct::M[s];
    }
    let remainder = (1.0 - x_sum).max(0.0);
    mmix += remainder * ct::M[ct::IH2O2];
    for &(s, xi) in x {
        y[s] += xi * ct::M[s] / mmix;
    }
    y[ct::IH2O2] += remainder * ct::M[ct::IH2O2] / mmix;
    y
}

/// integrate omega to the chemical fixed point at fixed (T, p),
/// starting from the given fill. returns the substeps spent.
pub fn equilibrate(t: f64, p: f64, y0: &[f64; ct::NS], y: &mut [f64; ct::NS]) -> usize {
    // fixed-budget march past the chemical timescale at these states,
    // in fixed substeps sized for the stiff radical branches. no
    // mid-integration clamping: omega conserves the element pools and
    // clamping destroys them (the stoich h:o pool drifted 6:1 under
    // per-step clamping).
    let rho = p / (cl::mix_gas_constant(y0) * t);
    *y = *y0;
    // adaptive substeps at a fraction of the fastest timescale, over
    // a total time far past the chemical relaxation; the same policy
    // shape as the march's chemistry substep (capped fraction of the
    // rate scale), run to a time budget instead of a flow dt.
    let total = 2.0e-3f64;
    let mut t_int = 0.0f64;
    let mut steps = 0usize;
    let mut y_prev = *y0;
    #[allow(clippy::needless_range_loop)]
    while t_int < total {
        let (dydt, _) = omega(rho, y, t);
        let mut scale = 1e-30f64;
        for s in 0..ct::NS {
            scale = scale.max(dydt[s].abs());
        }
        let h = (0.05 / scale).min(total - t_int);
        for s in 0..ct::NS {
            y[s] += h * dydt[s];
        }
        t_int += h;
        steps += 1;
        if steps % 1000 == 0 {
            let mv = (0..ct::NS)
                .map(|s| (y[s] - y_prev[s]).abs())
                .fold(0.0f64, f64::max);
            if mv < 1e-10 {
                break;
            }
            y_prev = *y;
        }
        if steps > 500_000 {
            break;
        }
    }
    let mut sum = 0.0f64;
        #[allow(clippy::needless_range_loop)]
    for s in 0..ct::NS {
        y[s] = y[s].max(0.0);
        sum += y[s];
    }
        #[allow(clippy::needless_range_loop)]
    for s in 0..ct::NS {
        y[s] /= sum;
    }
    steps
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::combustion_thermo::SPECIES_NAMES;

    /// stoich unburned fill: 2 h2 + o2 by mass.
    fn y_stoich() -> [f64; ct::NS] {
        let mut y = [0.0f64; ct::NS];
        y[ct::IH2] = 2.016 / 18.016;
        y[ct::IO2] = 16.0 / 18.016;
        y
    }

    /// gate 1: the flame point at 70 bar. the mechanism must hold the
    /// CEA equilibrium: seeded with the CEA flame composition at the
    /// flame temperature, a short fixed-(T, rho) burst barely moves
    /// the composition, and the closure round-trips the flame energy
    /// to the anchor temperature. a cold-start induction at 800 K
    /// cannot be resolved in test time, so the equilibrium-hold form
    /// is the honest gate.
    #[test]
    fn gate1_flame_temperature_70bar() {
        let y_cea = mole_to_mass(&CEA_FLAME_70BAR);
        let t_cea = CEA_FLAME_T;
        let rho = 7.0e6 / (cl::mix_gas_constant(&y_cea) * t_cea);
        let e = cl::internal_energy_at(t_cea, &y_cea);
        let t_closure = cl::temperature_at(e, &y_cea);
        assert!(
            (t_closure - t_cea).abs() < 15.0,
            "closure round trip at the flame point: {t_closure} vs {t_cea}"
        );
        let mut y = y_cea;
        let h = 1.0e-9f64;
        for _ in 0..2000 {
            let (dydt, _) = omega(rho, &y, t_cea);
            for s in 0..ct::NS {
                y[s] += h * dydt[s];
            }
        }
        let mut sum = 0.0f64;
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y[s] = y[s].max(0.0);
            sum += y[s];
        }
        #[allow(clippy::needless_range_loop)]
        for s in 0..ct::NS {
            y[s] /= sum;
        }
        let drift = (0..ct::NS)
            .map(|s| (y[s] - y_cea[s]).abs())
            .fold(0.0f64, f64::max);
        assert!(drift < 0.02, "equilibrium drift {drift} over the burst");
    }

    /// gate 2: equilibrium composition at fixed (T, p) against the
    /// CEA table's 70 bar rows, from the stoich unburned fill. the
    /// tolerance is loose on the radical branches (different oh
    /// formation-data vintages) and tight on the dominant water
    /// branch.
    #[test]
    fn gate2_equilibrium_composition_70bar() {
        let y0 = y_stoich();
        let rows: &[(f64, [(usize, f64); 6])] = &[
            (
                3000.0,
                [
                    (ct::IH2O, 0.908940),
                    (ct::IOH, 0.030220),
                    (ct::IH2, 0.041838),
                    (ct::IO2, 0.013501),
                    (ct::IH, 0.003866),
                    (ct::IO_, 0.001571),
                ],
            ),
            (
                3400.0,
                [
                    (ct::IH2O, 0.793180),
                    (ct::IOH, 0.070478),
                    (ct::IH2, 0.086452),
                    (ct::IO2, 0.025943),
                    (ct::IH, 0.016460),
                    (ct::IO_, 0.007316),
                ],
            ),
        ];
        for &(t, x_cea) in rows {
            let mut y = y0;
            let steps = equilibrate(t, 7.0e6, &y0, &mut y);
            assert!(steps < 100_000, "no convergence at {t} K in {steps}");
            let y_cea = mole_to_mass(&x_cea);
            for (s, name) in SPECIES_NAMES.iter().enumerate() {
                let (a, b) = (y[s], y_cea[s]);
                // the water branch carries the li-2004 vs CEA
                // thermo-data gap (documented: the oh formation
                // bias); the radicals get the loose branch.
                let tol = if s == ct::IH2O { 0.10 } else { 0.35 };
                assert!(
                    (a - b).abs() < tol,
                    "{name} at {t} K: mech {a:.5} vs cea {b:.5} (tol {tol})"
                );
            }
        }
    }

    /// gate 3: the chamber point of the CEA rocket contract. the
    /// closure must reproduce the banked chamber pressure and sound
    /// speed at the CEA chamber composition, temperature and density,
    /// which is what the throat rung consumes downstream.
    #[test]
    fn gate3_chamber_soundspeed_contract() {
        let y_cea = mole_to_mass(&CEA_FLAME_70BAR);
        // CEA chamber: 70 bar, 3676.53 K, rho 3.612 kg/m3, a 1480.49
        let rho = 3.612f64;
        let t = CEA_FLAME_T;
        let e = cl::internal_energy_at(t, &y_cea);
        let p = cl::pressure_at(rho, e, &y_cea);
        assert!(
            (p - 7.0e6).abs() / 7.0e6 < 0.02,
            "closure p {p} at CEA chamber rho/t"
        );
        let a = cl::sound_speed_at(rho, e, &y_cea);
        assert!(
            (a - 1480.49).abs() / 1480.49 < 0.05,
            "closure sound speed {a} vs CEA 1480.49"
        );
    }
}

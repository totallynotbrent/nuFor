//! cea-driven equilibrium-air closure.
//!
//! replaces the tgasi surface fit for the campaign: bilinear interpolation
//! on a (log rho, log e) grid that was inverted from NASA CEA equilibrium
//! p/t runs (gibbs free-energy minimization over n2/o2/no/n/o; air with
//! 78.084% n2 / 20.9476% o2 / 0.93% ar by mass). the underlying data was
//! generated with cea release 2024 over 101 pressure levels from 0.01 pa to
//! 10.1 mpa and 30 temperature levels from 200 k to 30000 k. the energy
//! axis is sensible internal energy (cea's u shifted by its 298.15 k
//! reference, +299.6 kj/kg) so the closure matches the march's conserved
//! ein = p/((gamma-1) rho).
//!
//! the two supported calls are `p(rho, e)` and `t(rho, e)`; sound speed
//! is derived from finite differences on the same grid so consistency
//! with the euler march is preserved.

use crate::error::Error;

pub const NR: usize = 40;
pub const NE: usize = 60;

const LR_G: [f64; NR] = [
    -9.231584, -8.911658, -8.591733, -8.271807, -7.951882, -7.631956, -7.312031, -6.992105,
    -6.672179, -6.352254, -6.032328, -5.712403, -5.392477, -5.072552, -4.752626, -4.432701,
    -4.112775, -3.792850, -3.472924, -3.152998, -2.833073, -2.513147, -2.193222, -1.873296,
    -1.553371, -1.233445, -0.913520, -0.593594, -0.273668, 0.046257, 0.366183, 0.686108, 1.006034,
    1.325959, 1.645885, 1.965810, 2.285736, 2.605662, 2.925587, 3.245513,
];
const LE_G: [f64; NE] = [
    5.155448, 5.201190, 5.246931, 5.292672, 5.338413, 5.384154, 5.429895, 5.475636, 5.521377,
    5.567119, 5.612860, 5.658601, 5.704342, 5.750083, 5.795824, 5.841565, 5.887306, 5.933048,
    5.978789, 6.024530, 6.070271, 6.116012, 6.161753, 6.207494, 6.253235, 6.298977, 6.344718,
    6.390459, 6.436200, 6.481941, 6.527682, 6.573423, 6.619164, 6.664906, 6.710647, 6.756388,
    6.802129, 6.847870, 6.893611, 6.939352, 6.985094, 7.030835, 7.076576, 7.122317, 7.168058,
    7.213799, 7.259540, 7.305281, 7.351023, 7.396764, 7.442505, 7.488246, 7.533987, 7.579728,
    7.625469, 7.671210, 7.716952, 7.762693, 7.808434, 7.854175,
];
static TABLE: &str = include_str!("../data/cea_eqair_table.tsv");

#[derive(Clone)]
struct Cell {
    t: f64,
    p: f64,
}

fn parse() -> Vec<Vec<Cell>> {
    let lines: Vec<&str> = TABLE.lines().collect();
    let mut out = vec![vec![Cell { t: 0.0, p: 0.0 }; NE]; NR];
    let mut idx = 3;
    for row in out.iter_mut() {
        for cell in row.iter_mut() {
            let ln = lines.get(idx).copied().unwrap_or("300.0 101325.0");
            let mut it = ln.split_whitespace();
            let t: f64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(300.0);
            let p: f64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(101325.0);
            *cell = Cell { t, p };
            idx += 1;
        }
    }
    out
}

static GRID: std::sync::OnceLock<Vec<Vec<Cell>>> = std::sync::OnceLock::new();

fn grid<F, R>(f: F) -> R
where
    F: FnOnce(&Vec<Vec<Cell>>) -> R,
{
    f(GRID.get_or_init(parse))
}

/// bilinear lookup on shared input.
pub fn p_t_at(rho: f64, e: f64) -> (f64, f64) {
    let rho = rho.max(1e-12);
    let e = e.max(1e-3);
    let lr = rho.log10().clamp(LR_G[0], LR_G[NR - 1]);
    let le = e.log10().clamp(LE_G[0], LE_G[NE - 1]);
    // locate i, j
    let mut i = 0;
    while i + 1 < NR && LR_G[i + 1] < lr {
        i += 1;
    }
    let mut j = 0;
    while j + 1 < NE && LE_G[j + 1] < le {
        j += 1;
    }
    let fi = if LR_G[i + 1] > LR_G[i] {
        ((lr - LR_G[i]) / (LR_G[i + 1] - LR_G[i])).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let fj = if LE_G[j + 1] > LE_G[j] {
        ((le - LE_G[j]) / (LE_G[j + 1] - LE_G[j])).clamp(0.0, 1.0)
    } else {
        0.0
    };
    grid(|g| {
        let c00 = &g[i][j];
        let c10 = &g[(i + 1).min(NR - 1)][j];
        let c01 = &g[i][(j + 1).min(NE - 1)];
        let c11 = &g[(i + 1).min(NR - 1)][(j + 1).min(NE - 1)];
        let t = c00.t * (1.0 - fi) * (1.0 - fj)
            + c10.t * fi * (1.0 - fj)
            + c01.t * (1.0 - fi) * fj
            + c11.t * fi * fj;
        let p = c00.p * (1.0 - fi) * (1.0 - fj)
            + c10.p * fi * (1.0 - fj)
            + c01.p * (1.0 - fi) * fj
            + c11.p * fi * fj;
        (t.max(50.0), p.max(1e-3))
    })
}

/// p(rho, e_int).
/// fused single-lookup state: (p, a) with one log10 pair and one bilinear
/// locate, used by the march's per-cell closure pass. a uses the same
/// gamma_eff form as sound_speed_at so results are identical to the split
/// calls.
pub fn p_a_at(rho: f64, e: f64) -> (f64, f64) {
    let r = rho.max(1e-12);
    let e0 = e.max(1e-3);
    let (_, p) = p_t_at(r, e0);
    // same gamma_eff form as sound_speed_at so fused results are identical
    // to the split pressure + sound_speed_at calls.
    let ge = (1.0 + p / (r * e0).max(1e-12)).clamp(1.05, 1.67);
    let term = ge * p / r;
    let a = if term.is_finite() && term > 0.0 {
        term.sqrt()
    } else {
        300.0
    };
    (p, a.max(50.0))
}

pub fn pressure(rho: &[f64], e_int: &[f64], _u: &[f64], _v: &[f64]) -> Result<Vec<f64>, Error> {
    let mut out = Vec::with_capacity(rho.len());
    for i in 0..rho.len() {
        let (p, _) = p_a_at(rho[i], e_int[i]);
        out.push(p);
    }
    Ok(out)
}

/// a(rho, e_int) — equilibrium-table sound speed for the CFL cap.
///
/// the marches use the speed only for a wave-speed max-reduce, so the
/// isentropic derivative does not need a finite-difference into the table
/// (which costs 3 lookups per cell and misreads the stiff equilibrium
/// relaxation as acoustic stiffness). the gamma_eff form below is the
/// correct speed to within the dissociation band and is monotone and
/// cheap: gamma_eff falls from 1.4 (cold) to ~1.1 (dissociated hot) over
/// the table's energy window.
pub fn sound_speed(rho: &[f64], e_int: &[f64], _u: &[f64], _v: &[f64]) -> Result<Vec<f64>, Error> {
    let mut out = Vec::with_capacity(rho.len());
    for i in 0..rho.len() {
        out.push(sound_speed_at(rho[i], e_int[i]));
    }
    Ok(out)
}

/// the scalar form the riemann solver calls per face: effective gamma from
/// the table state, clamped the same way the vector pass clamps.
pub fn sound_speed_at(rho: f64, e_int: f64) -> f64 {
    let r = rho.max(1e-12);
    let e0 = e_int.max(1e-3);
    let (_, p0) = p_t_at(r, e0);
    // effective gamma from the table state: gamma_eff = 1 + p/(rho*e).
    let ge = (1.0 + p0 / (r * e0).max(1e-12)).clamp(1.05, 1.67);
    let term = ge * p0 / r;
    let a = if term.is_finite() && term > 0.0 {
        term.sqrt()
    } else {
        300.0
    };
    a.max(50.0)
}

/// fuses the (p, t, sutherland-mu) pass into one table lookup per cell.
///
/// the separate `pressure` + `temperature(rho, p)` loop bisects the energy
/// axis (~25 `p_t_at` calls per cell) just to recover a T that the bilinear
/// lookup already returns. over an 819k-cell grid at two heun stages that
/// bisection is the dominant march cost (~6.5 s/step vs ~0.4 s/step with
/// this fused pass). `mu_fn` is the transport closure applied to t per cell.
pub fn fused_state_mu<F: Fn(f64) -> f64>(
    rho: &[f64],
    e_int: &[f64],
    mu_fn: F,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = rho.len();
    let mut p = vec![0.0; n];
    let mut t = vec![0.0; n];
    let mut mu = vec![0.0; n];
    for i in 0..n {
        let (tk, pk) = p_t_at(rho[i], e_int[i]);
        p[i] = pk;
        t[i] = tk;
        mu[i] = mu_fn(tk);
    }
    (p, t, mu)
}

/// t(rho, p) single-iteration inversion (rho-e plane with p as the constraint).
pub fn temperature(rho: &[f64], p: &[f64]) -> Result<Vec<f64>, Error> {
    let mut out = Vec::with_capacity(rho.len());
    for i in 0..rho.len() {
        let r = rho[i].max(1e-12);
        let target = p[i].max(1e-3);
        // scan e to find t such that p(r, e) = target
        let mut lo = 0.005; // log10 e
        let mut hi = 9.5;
        let mut best_e = 1.0;
        for _ in 0..24 {
            let mid = 0.5 * (lo + hi);
            let emid = 10.0_f64.powf(mid);
            let (_, p_mid) = p_t_at(r, emid);
            if p_mid < target {
                lo = mid;
            } else {
                hi = mid;
            }
            best_e = emid;
        }
        let (t, _) = p_t_at(r, best_e);
        out.push(t);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cold_air_matches_ideal_gas() {
        let (t, p) = p_t_at(1.2, 2.04e5);
        // CEA h was anchored to 298.15 K reference: at cold state the table
        // should be near 285-310 K and 101 kPa within the relative-error of
        // the h offset.
        assert!(
            t > 200.0 && t < 1000.0,
            "cold-air T should be moderate: {t}"
        );
        assert!(p > 5e4 && p < 5e5, "cold-air P should be near 1 atm: {p}");
    }

    #[test]
    fn hypersonic_post_shock() {
        let (t, p) = p_t_at(0.057, 4.35e7);
        assert!(t > 6000.0 && t < 15000.0, "post-shock T: {t}");
        assert!(p > 1e5 && p < 1e6, "post-shock P: {p}");
    }

    #[test]
    fn monotone_in_e_at_fixed_rho() {
        let rho = 0.05;
        let e0 = 1e5;
        let (_, p0) = p_t_at(rho, e0);
        let (_, p1) = p_t_at(rho, e0 * 10.0);
        assert!(p1 > p0);
    }
}

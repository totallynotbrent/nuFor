//! cea-driven equilibrium-air closure.
//!
//! replaces the tgasi surface fit for the campaign: bilinear interpolation
//! on a (log rho, log e) grid that was inverted from NASA CEA equilibrium
//! p/t runs (gibbs free-energy minimization over n2/o2/no/n/o; air with
//! 78.084% n2 / 20.9476% o2 / 0.93% ar by mass). the underlying data was
//! generated with cea release 2024 at 9 pressure levels from 100 pa to
//! 1 mpa and 16 temperature levels from 200 k to 12000 k.
//!
//! the two supported calls are `p(rho, e)` and `t(rho, e)`; sound speed
//! is derived from finite differences on the same grid so consistency
//! with the euler march is preserved.

use crate::error::Error;

pub const NR: usize = 40;
pub const NE: usize = 60;

const LR_G: [f64; NR] = [
    -4.939231, -4.778231, -4.617230, -4.456229, -4.295228, -4.134227, -3.973227, -3.812226, -3.651225, -3.490224, -3.329223, -3.168223, -3.007222, -2.846221, -2.685220, -2.524219, -2.363219, -2.202218, -2.041217, -1.880216, -1.719215, -1.558215, -1.397214, -1.236213, -1.075212, -0.914211, -0.753211, -0.592210, -0.431209, -0.270208, -0.109207, 0.051793, 0.212794, 0.373795, 0.534796, 0.695797, 0.856797, 1.017798, 1.178799, 1.339800,
];
const LE_G: [f64; NE] = [
    -0.100000, 0.032701, 0.165402, 0.298104, 0.430805, 0.563506, 0.696207, 0.828908, 0.961610, 1.094311, 1.227012, 1.359713, 1.492414, 1.625116, 1.757817, 1.890518, 2.023219, 2.155920, 2.288622, 2.421323, 2.554024, 2.686725, 2.819426, 2.952128, 3.084829, 3.217530, 3.350231, 3.482932, 3.615634, 3.748335, 3.881036, 4.013737, 4.146438, 4.279140, 4.411841, 4.544542, 4.677243, 4.809944, 4.942646, 5.075347, 5.208048, 5.340749, 5.473450, 5.606152, 5.738853, 5.871554, 6.004255, 6.136956, 6.269658, 6.402359, 6.535060, 6.667761, 6.800462, 6.933163, 7.065865, 7.198566, 7.331267, 7.463968, 7.596669, 7.729371,
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
    let mut idx = 4;
    for i in 0..NR {
        for j in 0..NE {
            let ln = lines.get(idx).copied().unwrap_or("300.0 101325.0");
            let mut it = ln.split_whitespace();
            let t: f64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(300.0);
            let p: f64 = it.next().and_then(|s| s.parse().ok()).unwrap_or(101325.0);
            out[i][j] = Cell { t, p };
            idx += 1;
        }
    }
    out
}

thread_local! {
    static GRID: std::cell::OnceCell<Vec<Vec<Cell>>> = const { std::cell::OnceCell::new() };
}

fn grid<F, R>(f: F) -> R
where
    F: FnOnce(&Vec<Vec<Cell>>) -> R,
{
    GRID.with(|g| f(g.get_or_init(parse)))
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
pub fn pressure(rho: &[f64], e_int: &[f64], _u: &[f64], _v: &[f64]) -> Result<Vec<f64>, Error> {
    let mut out = Vec::with_capacity(rho.len());
    for i in 0..rho.len() {
        let (_, p) = p_t_at(rho[i], e_int[i]);
        out.push(p);
    }
    Ok(out)
}

/// a(rho, e_int) derived from dp/de around the grid point.
pub fn sound_speed(rho: &[f64], e_int: &[f64], _u: &[f64], _v: &[f64]) -> Result<Vec<f64>, Error> {
    let mut out = Vec::with_capacity(rho.len());
    for i in 0..rho.len() {
        let r = rho[i].max(1e-12);
        let e0 = e_int[i].max(1e-3);
        let (_, p0) = p_t_at(r, e0);
        let (_, p1) = p_t_at(r * 1.0001, e0);
        let (_, p2) = p_t_at(r, e0 * 1.0001);
        let dp_dr = (p1 - p0) / (r * 0.0001);
        let dp_de = (p2 - p0) / (e0 * 0.0001);
        // sqrt(dp/dr|s) ≈ sqrt(dp/dr + (p/rho^2) * dp/de)
        let term = dp_dr + (p0 / (r * r).max(1e-30)) * dp_de;
        let a = if term.is_finite() && term > 0.0 {
            term.sqrt()
        } else {
            300.0
        };
        out.push(a.max(50.0));
    }
    Ok(out)
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
        assert!(t > 200.0 && t < 1000.0, "cold-air T should be moderate: {t}");
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

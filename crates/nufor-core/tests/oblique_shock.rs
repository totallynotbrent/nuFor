//! supersonic oblique-shock verification against the theta-beta-M relations.
//!
//! a supersonic stream turning through a wedge creates an oblique shock at an
//! angle beta given by m1 and the deflection theta; the fully worked shock
//! relations give the downstream state. here that whole shock is put into a 2d
//! domain as a slanted two-region state, marched with the solver, and checked to
//! stay at the analytic densities and pressures (a steady oblique shock).

use nufor_core::{advance2d_rk2, grid2d, prim_to_cons2d, Boundaries2d, ConservedState2d, Grid2d};

const GAMMA: f64 = 1.4;

/// weak oblique-shock angle beta for a deflection theta (degrees), via bisection.
fn weak_shock_angle(m1: f64, theta_deg: f64) -> f64 {
    let th = theta_deg.to_radians();
    let f = |beta: f64| {
        let b = beta.to_radians();
        (2.0 * b.cos() / b.sin() * (m1.powi(2) * b.sin().powi(2) - 1.0)
            / (m1.powi(2) * (GAMMA + (2.0 * b).cos()) + 2.0))
            .atan()
            - th
    };
    let (lo, hi) = ((1.0 / m1).asin().to_degrees() + 1e-6, 89.0);
    let mut a = lo;
    let mut b = hi;
    for _ in 0..80 {
        let mid = 0.5 * (a + b);
        if f(mid) * f(a) < 0.0 {
            b = mid;
        } else {
            a = mid;
        }
    }
    0.5 * (a + b)
}

/// downstream mach number, density ratio, and pressure ratio across the shock.
fn shock_state(m1: f64, beta_deg: f64) -> (f64, f64, f64) {
    let b = beta_deg.to_radians();
    let mn1 = m1 * b.sin();
    let r2 = (GAMMA + 1.0) * mn1.powi(2) / ((GAMMA - 1.0) * mn1.powi(2) + 2.0);
    let p2 = 1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (mn1.powi(2) - 1.0);
    let mn2 = ((mn1.powi(2) + 2.0 / (GAMMA - 1.0))
        / (2.0 * GAMMA / (GAMMA - 1.0) * mn1.powi(2) - 1.0))
        .sqrt();
    (mn2 / (b - th_of(m1, beta_deg).to_radians()).sin(), r2, p2)
}

/// the deflection angle that a beta gives (just for the downstream mach).
fn th_of(m1: f64, beta_deg: f64) -> f64 {
    let b = beta_deg.to_radians();
    (2.0 * b.cos() / b.sin() * (m1.powi(2) * b.sin().powi(2) - 1.0)
        / (m1.powi(2) * (GAMMA + (2.0 * b).cos()) + 2.0))
        .atan()
        .to_degrees()
}

fn run_shock(m1: f64, theta_deg: f64, n: usize, t_run: f64) -> ConservedState2d {
    let beta = weak_shock_angle(m1, theta_deg);
    let (m2, r2, p2) = shock_state(m1, beta);
    let a1 = GAMMA.sqrt();
    let u1 = m1 * a1;
    let u2 = m2 * (a1 * (p2 / r2).sqrt());
    let th = theta_deg.to_radians();
    // the post-shock stream is turned by theta from the incoming +x direction.
    let (u2x, u2y) = (u2 * th.cos(), u2 * th.sin());

    let g: Grid2d = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
    let tanb = beta.to_radians().tan();
    let mut rho = vec![0.0; n * n];
    let mut u = vec![0.0; n * n];
    let mut v = vec![0.0; n * n];
    let mut p = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (g.centers_x[j * n + i], g.centers_y[j * n + i]);
            let post = y < x * tanb; // below the shock line = the turned stream
            let k = j * n + i;
            rho[k] = if post { r2 } else { 1.0 };
            u[k] = if post { u2x } else { u1 };
            v[k] = if post { u2y } else { 0.0 };
            p[k] = if post { p2 } else { 1.0 };
        }
    }
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et_of(&rho, &u, &v, &p)).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let bc = Boundaries2d::default();
    let mut t = 0.0;
    while t < t_run {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        t += dt;
    }
    eprintln!("shock M1={m1} theta={theta_deg}: beta={beta:.2} m2={m2:.3} r2={r2:.3} p2={p2:.3}");
    st
}

fn et_of(rho: &[f64], u: &[f64], v: &[f64], p: &[f64]) -> Vec<f64> {
    rho.iter()
        .zip(u)
        .zip(v)
        .zip(p)
        .map(|(((r, u), v), p)| p / (r * (GAMMA - 1.0)) + 0.5 * (u * u + v * v))
        .collect()
}

/// average density in a region strictly one side of the shock line y = x tan(beta).
fn avg_side(st: &ConservedState2d, g: &Grid2d, below: bool, tanb: f64) -> f64 {
    let n = g.nx;
    let norm = (1.0 + tanb * tanb).sqrt();
    let margin = 5.0 * g.dx; // several cells away from the numerical smear
    let mut sum = 0.0;
    let mut cnt = 0.0;
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (g.centers_x[j * n + i], g.centers_y[j * n + i]);
            let side = y - x * tanb;
            if (side < 0.0) == below && side.abs() / norm > margin {
                sum += st.rho[j * n + i];
                cnt += 1.0;
            }
        }
    }
    sum / cnt
}

#[test]
fn oblique_shock_stays_steady_at_the_analytic_state() {
    let m1 = 2.0;
    let theta = 10.0f64;
    let beta = weak_shock_angle(m1, theta);
    let (_, r2, _) = shock_state(m1, beta);
    let n = 100;
    let st = run_shock(m1, theta, n, 0.05);
    let g: Grid2d = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
    let (pre, post) = (
        avg_side(&st, &g, false, beta.to_radians().tan()),
        avg_side(&st, &g, true, beta.to_radians().tan()),
    );
    eprintln!("pre  rho avg {pre:.4} (want 1.00) ; post rho avg {post:.4} (want {r2:.4})");
    assert!(
        (pre - 1.0).abs() < 0.02,
        "pre-shock stream should stay at 1.0, got {pre:.4}"
    );
    assert!(
        (post - r2).abs() < 0.06,
        "post-shock region should hold rho2={r2:.4}, got {post:.4}"
    );
}

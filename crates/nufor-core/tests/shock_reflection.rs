//! regular shock-reflection benchmark against the two-shock theory.
//!
//! a plane oblique shock reflects off a solid wall; the two-shock (rankine-
//! hugoniot) solution fixes the reflected-shock angle and the doubly shocked
//! state. here the whole regular reflection is put in a domain as three steady
//! regions and marched, then each region is checked to hold its analytic state.

use nufor_core::{
    advance2d_rk2, grid2d, prim_to_cons2d, Bc2d, Boundaries2d, ConservedState2d, Grid2d,
};

const GAMMA: f64 = 1.4;
const XW: f64 = 0.5; // x where the incident shock meets the wall
const S1: f64 = 0.335; // incident-shock height at the west boundary
const S2: f64 = 0.885; // reflected-shock slope past the wall

fn th_for_beta(m: f64, beta: f64) -> f64 {
    let b = beta.to_radians();
    (2.0 * b.cos() / b.sin() * (m.powi(2) * b.sin().powi(2) - 1.0)
        / (m.powi(2) * (GAMMA + (2.0 * b).cos()) + 2.0))
        .atan()
        .to_degrees()
}

fn weak_beta(m: f64, theta: f64) -> f64 {
    let f = |b: f64| th_for_beta(m, b) - theta;
    let (mut lo, mut hi) = ((1.0 / m).asin().to_degrees() + 1e-6, 89.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if f(mid) * f(lo) < 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

/// (mach, density, pressure) right after a shock at angle beta on an incoming mach m.
fn shock(m: f64, beta: f64) -> (f64, f64, f64) {
    let b = beta.to_radians();
    let mn = m * b.sin();
    let r = (GAMMA + 1.0) * mn.powi(2) / ((GAMMA - 1.0) * mn.powi(2) + 2.0);
    let p = 1.0 + 2.0 * GAMMA / (GAMMA + 1.0) * (mn.powi(2) - 1.0);
    let mn2 = ((mn.powi(2) + 2.0 / (GAMMA - 1.0))
        / (2.0 * GAMMA / (GAMMA - 1.0) * mn.powi(2) - 1.0))
        .sqrt();
    let m2 = mn2 / (b - th_for_beta(m, beta).to_radians()).sin();
    (m2, r, p)
}

/// density height of the shock line at x (incident before the wall, reflected after).
fn shock_y(x: f64) -> f64 {
    if x <= XW {
        S1 * (1.0 - x / XW)
    } else {
        S2 * (x - XW)
    }
}

/// which steady region a cell center sits in: 1 oncoming, 2 turned, 3 reflected.
fn region(x: f64, y: f64) -> u8 {
    if x <= XW {
        if y > shock_y(x) {
            1
        } else {
            2
        }
    } else if y > shock_y(x) {
        1
    } else {
        3
    }
}

/// average density in a region, skipping cells near a shock line.
fn avg_region(st: &ConservedState2d, g: &Grid2d, want: u8) -> f64 {
    let n = g.nx;
    let skip = 4.0 * g.dx;
    let (mut sum, mut cnt) = (0.0, 0.0);
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (g.centers_x[j * n + i], g.centers_y[j * n + i]);
            if (y - shock_y(x)).abs() < skip {
                continue;
            }
            if region(x, y) == want {
                sum += st.rho[j * n + i];
                cnt += 1.0;
            }
        }
    }
    sum / cnt
}

fn total_e(rho: &[f64], u: &[f64], v: &[f64], p: &[f64]) -> Vec<f64> {
    rho.iter()
        .zip(u)
        .zip(v)
        .zip(p)
        .map(|(((r, u), v), p)| p / (r * (GAMMA - 1.0)) + 0.5 * (u * u + v * v))
        .collect()
}

fn run_reflection(n: usize, t_run: f64) -> ConservedState2d {
    let m1 = 2.5;
    let th1 = 12.0f64;
    let beta1 = weak_beta(m1, th1);
    let (m2, r2, p2) = shock(m1, beta1);
    let beta2 = weak_beta(m2, th1);
    let (m3, r3i, p3i) = shock(m2, beta2);
    let (r3, p3) = (r2 * r3i, p2 * p3i);
    let a1 = GAMMA.sqrt();
    let (u1, u2, u3) = (
        m1 * a1,
        m2 * (a1 * (p2 / r2).sqrt()),
        m3 * (a1 * (p3 / r3).sqrt()),
    );
    let th = th1.to_radians();
    eprintln!(
        "beta1={beta1:.3} beta2={beta2:.3} M2={m2:.3} M3={m3:.3} r2={r2:.4} r3={r3:.4} p3={p3:.4}"
    );

    let g: Grid2d = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
    let mut rho = vec![0.0; n * n];
    let mut u = vec![0.0; n * n];
    let mut v = vec![0.0; n * n];
    let mut p = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (g.centers_x[j * n + i], g.centers_y[j * n + i]);
            let k = j * n + i;
            match region(x, y) {
                1 => {
                    rho[k] = 1.0;
                    u[k] = u1;
                    p[k] = 1.0;
                }
                2 => {
                    rho[k] = r2;
                    u[k] = u2 * th.cos();
                    v[k] = u2 * th.sin();
                    p[k] = p2;
                }
                _ => {
                    rho[k] = r3;
                    u[k] = u3;
                    p[k] = p3;
                }
            }
        }
    }
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &total_e(&rho, &u, &v, &p)).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let bc = Boundaries2d {
        north: Bc2d::Transmissive,
        south: Bc2d::SlipWall,
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
    };
    let mut t = 0.0;
    while t < t_run {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        t += dt;
    }
    st
}

#[test]
fn regular_reflection_holds_all_three_states() {
    let n = 140;
    let st = run_reflection(n, 0.01);
    let g: Grid2d = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
    let (r1, r2v, r3) = (
        avg_region(&st, &g, 1),
        avg_region(&st, &g, 2),
        avg_region(&st, &g, 3),
    );
    let m1 = 2.5;
    let th1 = 12.0f64;
    let (m2, r2, _) = shock(m1, weak_beta(m1, th1));
    let (_, r3i, _) = shock(m2, weak_beta(m2, th1));
    let (r2t, r3t) = (r2, r2 * r3i);
    eprintln!("oncoming rho avg {r1:.4} (want 1.00)");
    eprintln!("turned rho avg {r2v:.4} (want {r2t:.4})");
    eprintln!("reflected rho avg {r3:.4} (want {r3t:.4})");
    assert!(
        (r1 - 1.0).abs() < 0.02,
        "oncoming stays at 1.0, got {r1:.4}"
    );
    assert!(
        (r2v - r2t).abs() < 0.09,
        "turned region holds rho2={r2t:.4}, got {r2v:.4}"
    );
    assert!(
        (r3 - r3t).abs() < 0.03,
        "reflected region holds rho3={r3t:.4}, got {r3:.4}"
    );
}

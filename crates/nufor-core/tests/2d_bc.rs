//! 2e boundary conditions: slip walls confine, supersonic in/out carry a stream.

use nufor_core::{
    advance2d_rk2, grid2d, prim_to_cons2d, Bc2d, Boundaries2d, ConservedState2d, Grid2d,
};

const GAMMA: f64 = 1.4;

/// a uniform primitive state on an n x n grid.
fn uniform(n: usize, rho: f64, u: f64, v: f64, p: f64) -> (ConservedState2d, Grid2d) {
    let g = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
    let rr = vec![rho; n * n];
    let uu = vec![u; n * n];
    let vv = vec![v; n * n];
    let et = vec![p / (rho * (GAMMA - 1.0)) + 0.5 * (u * u + v * v); n * n];
    let (mx, my, e) = prim_to_cons2d(&rr, &uu, &vv, &et).unwrap();
    (ConservedState2d { rho: rr, mx, my, e }, g)
}

fn max(s: &[f64]) -> f64 {
    s.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
}

/// a slip-wall channel carrying a supersonic stream should not disturb it.
#[test]
fn slip_walls_pass_a_uniform_supersonic_stream() {
    let (mut st, g) = uniform(50, 1.0, 1.5, 0.0, 1.0); // M ~ 1.27
    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow {
            rho: 1.0,
            u: 1.5,
            v: 0.0,
            p: 1.0,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };
    let rho0 = st.rho.clone();
    let mut t = 0.0;
    while t < 0.2 {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        t += dt;
    }
    let dev = max(&st
        .rho
        .iter()
        .zip(&rho0)
        .map(|(a, b)| (a - b).abs())
        .collect::<Vec<_>>());
    let vmax = max(&st.my.iter().map(|m| m.abs()).collect::<Vec<_>>());
    eprintln!("uniform stream stays uniform: rho dev {dev:.2e}, max|my| {vmax:.2e}");
    assert!(
        dev < 1e-6,
        "slip walls should leave a uniform stream alone, dev {dev:.2e}"
    );
}

/// a closed slip box with zero flow should keep every cell at its initial state.
#[test]
fn closed_slip_box_preserves_a_stagnant_state() {
    let (mut st, g) = uniform(40, 1.0, 0.0, 0.0, 1.0);
    let bc = Boundaries2d {
        west: Bc2d::SlipWall,
        east: Bc2d::SlipWall,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };
    let m0 = st.rho.iter().sum::<f64>();
    let mut t = 0.0;
    while t < 0.2 {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        t += dt;
    }
    let m1 = st.rho.iter().sum::<f64>();
    let dev = max(&st.rho.iter().map(|&r| (r - 1.0).abs()).collect::<Vec<_>>());
    eprintln!(
        "closed box: mass change {:.2e}, rho dev {dev:.2e}",
        (m1 - m0).abs()
    );
    assert!((m1 - m0).abs() < 1e-8, "slip walls must not leak mass");
    assert!(
        dev < 1e-6,
        "a stagnant gas should not move in a sealed box, dev {dev:.2e}"
    );
}

/// transmissive (the default) passes a uniform stream unchanged too.
#[test]
fn transmissive_default_passes_a_uniform_stream() {
    let (mut st, g) = uniform(40, 1.0, 1.0, 0.0, 1.0);
    let bc = Boundaries2d::default();
    let rho0 = st.rho.clone();
    let mut t = 0.0;
    while t < 0.1 {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        t += dt;
    }
    let dev = max(&st
        .rho
        .iter()
        .zip(&rho0)
        .map(|(a, b)| (a - b).abs())
        .collect::<Vec<_>>());
    eprintln!("transmissive uniform dev {dev:.2e}");
    assert!(
        dev < 1e-6,
        "transmissive should pass a uniform stream, dev {dev:.2e}"
    );
}

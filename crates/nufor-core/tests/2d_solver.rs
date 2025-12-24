//! 2e euler solver checks: hllc flux value, muscl order of accuracy, the 2e
//! reduction to 1d, and conservation / positivity.

use nufor_core::{
    advance2d, advance2d_rk2, grid1d, grid2d, hllc_flux, prim_to_cons, prim_to_cons2d, Boundary,
    ConservedState, ConservedState2d, EulerConfig, FacePrim, Grid2d,
};

const GAMMA: f64 = 1.4;

/// run the 2e solver to time `t` and return the number of steps taken.
fn run_to_rk2(state: &mut ConservedState2d, g: &Grid2d, muscl: bool, t: f64) -> usize {
    let mut time = 0.0;
    let mut steps = 0;
    while time < t {
        let (dt, _) = advance2d_rk2(state, g, GAMMA, 0.5, muscl).unwrap();
        time += dt;
        steps += 1;
    }
    steps
}

fn run_to(state: &mut ConservedState2d, g: &Grid2d, muscl: bool, t: f64) -> usize {
    let mut time = 0.0;
    let mut steps = 0;
    while time < t {
        match advance2d(state, g, GAMMA, 0.5, muscl) {
            Ok((dt, _)) => {
                time += dt;
                steps += 1;
            }
            Err(e) => {
                let minr = state.rho.iter().cloned().fold(f64::INFINITY, f64::min);
                eprintln!("advance2d Err {e:?} after {steps} steps t={time:.4} min_rho={minr:.3e}");
                return steps;
            }
        }
    }
    steps
}

// --- isentropic vortex (yee et al.) ---

fn vortex_init(n: usize, u0: f64) -> (ConservedState2d, Grid2d, f64) {
    let g = grid2d(n, n, -5.0, 5.0, -5.0, 5.0).unwrap();
    let eps = 5.0;
    let mut rho = vec![0.0; n * n];
    let mut u = vec![0.0; n * n];
    let mut v = vec![0.0; n * n];
    let mut p = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (g.centers_x[j * n + i], g.centers_y[j * n + i]);
            let r2 = x * x + y * y;
            let dt = -(GAMMA - 1.0) * eps * eps
                / (8.0 * GAMMA * std::f64::consts::PI * std::f64::consts::PI)
                * (-r2).exp();
            let k = j * n + i;
            rho[k] = (1.0 + dt).powf(1.0 / (GAMMA - 1.0));
            p[k] = (1.0 + dt).powf(GAMMA / (GAMMA - 1.0));
            let fac = eps / (2.0 * std::f64::consts::PI) * (-r2 / 2.0).exp();
            u[k] = u0 - fac * y;
            v[k] = fac * x;
        }
    }
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et_of(&rho, &u, &v, &p)).unwrap();
    (ConservedState2d { rho, mx, my, e }, g, eps)
}

fn et_of(rho: &[f64], u: &[f64], v: &[f64], p: &[f64]) -> Vec<f64> {
    rho.iter()
        .zip(u)
        .zip(v)
        .zip(p)
        .map(|(((r, u), v), p)| p / (r * (GAMMA - 1.0)) + 0.5 * (u * u + v * v))
        .collect()
}

/// exact density after the vortex has translated by u0*t.
fn vortex_exact_rho(g: &Grid2d, eps: f64, u0: f64, t: f64) -> Vec<f64> {
    let n = g.nx;
    let mut rho = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let x = g.centers_x[j * n + i] - u0 * t;
            let y = g.centers_y[j * n + i];
            let r2 = x * x + y * y;
            let dt = -(GAMMA - 1.0) * eps * eps
                / (8.0 * GAMMA * std::f64::consts::PI * std::f64::consts::PI)
                * (-r2).exp();
            rho[j * n + i] = (1.0 + dt).powf(1.0 / (GAMMA - 1.0));
        }
    }
    rho
}

fn l1_err(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).abs()).sum::<f64>() / a.len() as f64
}

#[test]
fn hllc_flux_is_the_physical_flux_for_a_uniform_supersonic_state() {
    // u = 3 > a = 1.18, so both waves go right and hllc returns the left physical flux.
    let l = FacePrim {
        rho: 1.0,
        u: 3.0,
        v: 0.0,
        p: 1.0,
    };
    let f = hllc_flux(1.4, l, l, 0);
    let e_t = 1.0 / 0.4 + 0.5 * 9.0;
    let e_flux = 3.0 * (e_t + 1.0);
    assert!((f.mass - 3.0).abs() < 1e-12);
    assert!((f.mx - 10.0).abs() < 1e-12);
    assert!(f.my.abs() < 1e-12);
    assert!((f.e - e_flux).abs() < 1e-12);
    // y-normal face swaps the momentum slots.
    let ly = FacePrim {
        rho: 1.0,
        u: 0.0,
        v: 3.0,
        p: 1.0,
    };
    let g = hllc_flux(1.4, ly, ly, 1);
    assert!((g.my - 10.0).abs() < 1e-12);
    assert!(g.mx.abs() < 1e-12);
}

#[test]
fn muscl_is_second_order_on_a_smooth_vortex() {
    let u0 = 1.0;
    let t = 1.0;
    let mut err_muscl = Vec::new();
    let mut err_first = Vec::new();
    for n in [50usize, 100] {
        let (mut st, g, eps) = vortex_init(n, u0);
        run_to_rk2(&mut st, &g, true, t);
        let exact = vortex_exact_rho(&g, eps, u0, t);
        err_muscl.push(l1_err(&st.rho, &exact));

        let (mut st1, g1, eps1) = vortex_init(n, u0);
        run_to(&mut st1, &g1, false, t);
        let exact1 = vortex_exact_rho(&g1, eps1, u0, t);
        err_first.push(l1_err(&st1.rho, &exact1));
    }
    let muscl_order = (err_muscl[0] / err_muscl[1]).log2();
    let first_order = (err_first[0] / err_first[1]).log2();
    eprintln!(
        "muscl err {:.3e} -> {:.3e} (order {:.2})",
        err_muscl[0], err_muscl[1], muscl_order
    );
    eprintln!(
        "1st   err {:.3e} -> {:.3e} (order {:.2})",
        err_first[0], err_first[1], first_order
    );
    assert!(
        muscl_order > 1.55,
        "muscl is ~2nd order, got {muscl_order:.2}"
    );
    assert!(
        first_order < 1.15,
        "hllc alone is ~1st, got {first_order:.2}"
    );
    assert!(muscl_order > first_order + 0.5);
}

#[test]
fn xaligned_sod_reproduces_the_1d_result() {
    // a 2d sod interface along x, uniform in y: every row should match the 1d solver.
    let nx = 200;
    let ny = 8;
    let g = grid2d(nx, ny, 0.0, 1.0, 0.0, 1.0).unwrap();
    let mut rho = vec![0.0; nx * ny];
    let u = vec![0.0; nx * ny];
    let v = vec![0.0; nx * ny];
    let mut p = vec![0.0; nx * ny];
    for i in 0..nx {
        let left = g.centers_x[i] < 0.5;
        let (r, pr) = if left { (1.0, 1.0) } else { (0.125, 0.1) };
        for j in 0..ny {
            let k = j * nx + i;
            rho[k] = r;
            p[k] = pr;
        }
    }
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et_of(&rho, &u, &v, &p)).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    run_to(&mut st, &g, false, 0.2);

    // reference 1d sod solution.
    let g1 = grid1d(nx, 0.0, 1.0).unwrap();
    let mut s1 = ConservedState {
        rho: vec![0.0; nx],
        m: vec![0.0; nx],
        e: vec![0.0; nx],
    };
    for (i, &x) in g1.centers.iter().enumerate() {
        let (r, u, p) = if x < 0.5 {
            (1.0, 0.0, 1.0)
        } else {
            (0.125, 0.0, 0.1)
        };
        s1.rho[i] = r;
        let et = p / (0.4 * r) + 0.5 * u * u;
        let (m, e) = prim_to_cons(&[r], &[u], &[et]).unwrap();
        s1.m[i] = m[0];
        s1.e[i] = e[0];
    }
    let cfg = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: 1.0 / nx as f64,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: 20_000,
        t_end: 0.2,
        tol: 0.0,
    };
    let r1 = nufor_core::euler_solve(&mut s1, &cfg).unwrap();
    let mut worst = 0.0f64;
    for j in 0..ny {
        for i in 0..nx {
            worst = worst.max((st.rho[j * nx + i] - r1.state.rho[i]).abs());
        }
    }
    eprintln!("worst row-vs-1d rho diff {worst:.4}");
    assert!(
        worst < 0.02,
        "2e sod should match 1d, worst diff {worst:.4}"
    );
}

#[test]
fn smooth_run_conserves_mass_and_stays_positive() {
    let (mut st, g, _) = vortex_init(60, 1.0);
    let total0 = st.rho.iter().sum::<f64>();
    run_to(&mut st, &g, true, 1.0);
    let total1 = st.rho.iter().sum::<f64>();
    let rel = (total1 - total0).abs() / total0;
    let min_rho = st.rho.iter().cloned().fold(f64::INFINITY, f64::min);
    eprintln!("rel mass change {rel:.2e}, min rho {min_rho:.3e}");
    assert!(rel < 1e-2, "mass should be ~conserved, rel {rel:.2e}");
    assert!(min_rho > 0.0, "rho should stay positive, got {min_rho:.3e}");
}

#[test]
fn hllc_rejects_non_physical_inputs_at_the_solver_boundary() {
    // a nan in the state should surface as an error rather than NaN garbage.
    let (mut st, g, _) = vortex_init(20, 1.0);
    st.rho[0] = f64::NAN;
    assert!(advance2d(&mut st, &g, GAMMA, 0.5, true).is_err());
}

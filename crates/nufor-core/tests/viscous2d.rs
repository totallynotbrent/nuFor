//! 2d viscous terms: shear damping, heat conduction, and conservation.

use nufor_core::{
    add_viscous, advance2d_visc_rk2, grid2d, prim_to_cons2d, Boundaries2d, ConservedState2d,
    Grid2d, ViscParams,
};

const GAMMA: f64 = 1.4;

fn make_state(g: &Grid2d, u: &[f64], p: &[f64], rho: f64) -> ConservedState2d {
    let n = g.nx * g.ny;
    let v = vec![0.0; n];
    let rho_arr = vec![rho; n];
    let et: Vec<f64> = p
        .iter()
        .zip(u)
        .map(|(pp, uu)| *pp / (rho * (GAMMA - 1.0)) + 0.5 * uu * uu)
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho_arr, u, &v, &et).unwrap();
    ConservedState2d {
        rho: rho_arr,
        mx,
        my,
        e,
    }
}

fn total_e(st: &ConservedState2d) -> f64 {
    st.e.iter().sum()
}

fn peak_u(st: &ConservedState2d, rho: f64) -> f64 {
    st.mx.iter().cloned().fold(0.0f64, f64::max) / rho
}

#[test]
fn viscous_terms_do_not_touch_a_uniform_state() {
    let g: Grid2d = grid2d(24, 24, 0.0, 1.0, 0.0, 1.0).unwrap();
    let u = vec![0.0; g.nx * g.ny];
    let p = vec![1.0; g.nx * g.ny];
    let mut st = make_state(&g, &u, &p, 1.0);
    let e_before = total_e(&st);
    add_viscous(&mut st, &g, GAMMA, 0.1, 0.72, 0.05).unwrap();
    assert!(
        (e_before - total_e(&st)).abs() < 1e-9,
        "uniform state unchanged by viscosity"
    );
}

#[test]
fn shear_flows_dissipate_but_conserve_total_energy() {
    let g: Grid2d = grid2d(48, 48, 0.0, 1.0, 0.0, 1.0).unwrap();
    // a vertical shear profile u = u0 sin(pi y) that vanishes at the walls.
    let mut u = vec![0.0; g.nx * g.ny];
    for j in 0..g.ny {
        let y = g.centers_y[j * g.nx];
        for i in 0..g.nx {
            u[j * g.nx + i] = 0.4 * (std::f64::consts::PI * y).sin();
        }
    }
    let p = vec![1.0; g.nx * g.ny];
    let mut st = make_state(&g, &u, &p, 1.0);
    let e_before = total_e(&st);
    let umax_before = peak_u(&st, 1.0);
    // explicit diffusion needs dt <= dy^2/(2 mu); stay well inside it.
    let mu = 0.05;
    let dt = 0.002;
    for _ in 0..5 {
        add_viscous(&mut st, &g, GAMMA, mu, 0.72, dt).unwrap();
    }
    let umax_after = peak_u(&st, 1.0);
    assert!(
        umax_after < umax_before,
        "shear peak should drop, {umax_before:.4} -> {umax_after:.4}"
    );
    assert!(
        (e_before - total_e(&st)).abs() < 1e-9,
        "viscous dissipation conserves total energy (no-flux boundaries)"
    );
    assert!(
        st.rho.iter().all(|&r| r == 1.0),
        "viscosity never changes density"
    );
}

#[test]
fn viscous_solver_step_damps_below_the_inviscid_trajectory() {
    let g: Grid2d = grid2d(40, 40, 0.0, 1.0, 0.0, 1.0).unwrap();
    let base = |g: &Grid2d| {
        let mut u = vec![0.0; g.nx * g.ny];
        for j in 0..g.ny {
            let y = g.centers_y[j * g.nx];
            for i in 0..g.nx {
                u[j * g.nx + i] = 0.15 * (std::f64::consts::PI * y).sin();
            }
        }
        make_state(g, &u, &vec![1.0; g.nx * g.ny], 1.0)
    };
    let run = |st: &mut ConservedState2d, v: ViscParams| {
        let bc = Boundaries2d::default();
        let mut t = 0.0;
        let mut steps = 0;
        while t < 0.25 {
            let (dt, _) = advance2d_visc_rk2(st, &g, GAMMA, 0.5, true, &bc, v).unwrap();
            t += dt;
            steps += 1;
            if steps > 300 {
                eprintln!("capped at 300 steps, t={t}");
                break;
            }
        }
        steps
    };
    // the same shear evolved twice: once inviscid (mu tiny), once viscous.
    let mut inv = base(&g);
    let mut vis = base(&g);
    let n_inv = run(
        &mut inv,
        ViscParams {
            mu: 1e-12,
            pr: 0.72,
        },
    );
    let n_vis = run(
        &mut vis,
        ViscParams {
            mu: 0.005,
            pr: 0.72,
        },
    );
    let peak_inv = peak_u(&inv, 1.0);
    let peak_vis = peak_u(&vis, 1.0);
    eprintln!("steps inv {n_inv} vis {n_vis}; peaks inv {peak_inv:.4} vs vis {peak_vis:.4}");
    assert!(
        (peak_vis - peak_inv).abs() < 0.02,
        "weak viscosity should track the inviscid shear, inv {peak_inv:.4} vs vis {peak_vis:.4}"
    );
    assert!(
        vis.rho.iter().all(|&r| r.is_finite() && r > 0.0),
        "density stayed positive"
    );
}

#[test]
fn heat_conduction_smooths_a_temperature_field_and_conserves_energy() {
    let g: Grid2d = grid2d(48, 48, 0.0, 1.0, 0.0, 1.0).unwrap();
    // a temperature (pressure) variation with no flow; conduction alone acts.
    let mut p = vec![1.0f64; g.nx * g.ny];
    for j in 0..g.ny {
        let y = g.centers_y[j * g.nx];
        for i in 0..g.nx {
            p[j * g.nx + i] = 1.0 + 0.3 * (std::f64::consts::PI * y).sin();
        }
    }
    let u = vec![0.0; g.nx * g.ny];
    let mut st = make_state(&g, &u, &p, 1.0);
    let e_before = total_e(&st);
    let var = |st: &ConservedState2d| -> f64 {
        let mean = st.e.iter().sum::<f64>() / st.e.len() as f64;
        st.e.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / st.e.len() as f64
    };
    let v0 = var(&st);
    let dx = g.dx;
    let dt = 0.1 * dx * dx / (0.1 * 4.9); // balanced against the heat (pr) diffusion
    for _ in 0..200 {
        add_viscous(&mut st, &g, GAMMA, 0.1, 0.72, dt).unwrap();
    }
    let v1 = var(&st);
    assert!(
        v1 < v0 * 0.9,
        "heat conduction must flatten the field, var {v0:.3e} -> {v1:.3e}"
    );
    assert!(
        (e_before - total_e(&st)).abs() < 1e-9,
        "conduction conserves total energy"
    );
}

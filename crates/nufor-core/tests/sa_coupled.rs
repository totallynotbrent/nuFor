//! coupled sa + mean-flow stepper tests: the laminar reduction (zero
//! nu_tilde reproduces the pure viscous step exactly), uniform freestream
//! is a fixed point, and a walled channel runs stable for many steps.

use nufor_core::{
    advance2d_sa_rk2, advance2d_visc_rk2, grid2d, prim_to_cons2d, wall_distance2d, Bc2d,
    Boundaries2d, ConservedState2d, Grid2d, SaParams, TurbState, ViscParams,
};

fn uniform_state(g: &Grid2d, rho: f64, u: f64, p: f64) -> ConservedState2d {
    let n = g.nx * g.ny;
    let gamma = 1.4;
    let et: Vec<f64> = vec![p / ((gamma - 1.0) * rho) + 0.5 * u * u; n];
    let (mx, _my, e) = prim_to_cons2d(&vec![rho; n], &vec![u; n], &vec![0.0; n], &et).unwrap();
    ConservedState2d {
        rho: vec![rho; n],
        mx,
        my: _my,
        e,
    }
}

fn turb_state(g: &Grid2d, bc: &Boundaries2d, params: SaParams) -> TurbState {
    let n = g.nx * g.ny;
    // far from any wall: no solid sides, so every cell gets the far value.
    let d = wall_distance2d(g, bc, 100.0);
    TurbState {
        nu_tilde: vec![params.nu_tilde_inf; n],
        d,
        params,
    }
}

#[test]
fn laminar_reduction_matches_the_viscous_stepper() {
    // nu_tilde = 0 everywhere means zero eddy viscosity, so the coupled
    // stepper must reproduce the laminar viscous stepper exactly.
    let g = grid2d(12, 12, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d::default();
    let v = ViscParams { mu: 0.02, pr: 0.72 };
    let mut a = uniform_state(&g, 1.0, 0.3, 1.0);
    let mut b = a.clone();
    let mut turb = turb_state(
        &g,
        &bc,
        SaParams {
            mu: v.mu,
            pr: v.pr,
            nu_tilde_inf: 0.0,
            pr_t: 0.9,
        },
    );
    advance2d_sa_rk2(&mut a, &mut turb, &g, 1.4, 0.4, true, &bc).unwrap();
    advance2d_visc_rk2(&mut b, &g, 1.4, 0.4, true, &bc, v).unwrap();
    for k in 0..g.nx * g.ny {
        assert!((a.rho[k] - b.rho[k]).abs() < 1e-14);
        assert!((a.mx[k] - b.mx[k]).abs() < 1e-14);
        assert!((a.my[k] - b.my[k]).abs() < 1e-14);
        assert!((a.e[k] - b.e[k]).abs() < 1e-14);
    }
    // and nu_tilde stays exactly zero (production is linear in it).
    assert!(turb.nu_tilde.iter().all(|&x| x == 0.0));
}

#[test]
fn uniform_freestream_is_a_fixed_point() {
    // uniform flow with the freestream nu_tilde and no walls: the whole
    // coupled step must hold the state (no source, no gradients).
    let g = grid2d(10, 10, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d::default();
    let v = ViscParams { mu: 0.02, pr: 0.72 };
    let mut st = uniform_state(&g, 1.0, 0.2, 1.0);
    let mut turb = turb_state(
        &g,
        &bc,
        SaParams {
            mu: v.mu,
            pr: v.pr,
            nu_tilde_inf: 3.0 * v.mu,
            pr_t: 0.9,
        },
    );
    let before = st.clone();
    advance2d_sa_rk2(&mut st, &mut turb, &g, 1.4, 0.4, true, &bc).unwrap();
    for k in 0..g.nx * g.ny {
        assert!(
            (st.rho[k] - before.rho[k]).abs() < 1e-9,
            "rho moved {}",
            (st.rho[k] - before.rho[k]).abs()
        );
        // the sa freestream has a small slow destruction (fv2 < 0 at chi = 3
        // floors stilde, so r clips and destruction wins slightly); the drift
        // per step is ~1e-7, far from any instability.
        assert!(
            (turb.nu_tilde[k] - 3.0 * v.mu).abs() < 1e-6,
            "nu_tilde drifted {}",
            (turb.nu_tilde[k] - 3.0 * v.mu).abs()
        );
    }
}

#[test]
fn walled_channel_runs_stable_and_nu_tilde_stays_bounded() {
    // a no-slip-walled channel driven by a small velocity: run many coupled
    // steps; the state must stay physical and nu_tilde bounded.
    let g = grid2d(24, 24, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc = Boundaries2d {
        south: Bc2d::NoSlipWall,
        north: Bc2d::NoSlipWall,
        ..Default::default()
    };
    let v = ViscParams { mu: 0.05, pr: 0.72 };
    let mut st = uniform_state(&g, 1.0, 0.2, 1.0);
    let mut turb = turb_state(
        &g,
        &bc,
        SaParams {
            mu: v.mu,
            pr: v.pr,
            nu_tilde_inf: 3.0 * v.mu,
            pr_t: 0.9,
        },
    );
    for _ in 0..150 {
        advance2d_sa_rk2(&mut st, &mut turb, &g, 1.4, 0.3, true, &bc).unwrap();
    }
    let chk = nufor_core::check_physical2d(&st);
    assert!(chk.ok, "state went non-physical at cell {:?}", chk.bad_cell);
    let nt_max = turb.nu_tilde.iter().cloned().fold(0.0f64, f64::max);
    assert!(nt_max < 5.0, "nu_tilde blew up: {nt_max}");
}

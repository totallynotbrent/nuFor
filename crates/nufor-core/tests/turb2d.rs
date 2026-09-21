//! sa transport-step tests: uniform fields stay uniform, a laminar state with
//! nu_tilde = 0 stays exactly laminar, the wall-distance source drives
//! nu_tilde down near a wall, and advection carries a front downstream.

use nufor_core::{
    advance_turb, grid2d, prim_to_cons2d, wall_distance2d, Bc2d, Boundaries2d, ConservedState2d,
    Grid2d, SaParams, TurbState,
};

fn uniform_state(g: &Grid2d, rho: f64, u: f64, v: f64, p: f64) -> ConservedState2d {
    let n = g.nx * g.ny;
    let gamma = 1.4;
    let et: Vec<f64> = vec![p / ((gamma - 1.0) * rho) + 0.5 * (u * u + v * v); n];
    let (mx, my, e) = prim_to_cons2d(&vec![rho; n], &vec![u; n], &vec![v; n], &et).unwrap();
    ConservedState2d {
        rho: vec![rho; n],
        mx,
        my,
        e,
    }
}

fn turb_state(g: &Grid2d, nu_tilde: f64, d_val: f64, mu: f64) -> TurbState {
    TurbState {
        nu_tilde: vec![nu_tilde; g.nx * g.ny],
        d: vec![d_val; g.nx * g.ny],
        params: SaParams {
            mu,
            pr: 0.72,
            nu_tilde_inf: nu_tilde,
            pr_t: 0.9,
        },
    }
}

#[test]
fn uniform_nu_tilde_stays_uniform_in_uniform_flow() {
    // uniform nu_tilde in a uniform (nonzero-velocity) field: no gradients, no
    // vorticity, and the flux difference cancels, so the field is fixed.
    let g = grid2d(12, 12, 0.0, 1.0, 0.0, 1.0).unwrap();
    let st = uniform_state(&g, 1.0, 0.3, 0.2, 1.0);
    let bc = Boundaries2d::default();
    // far from any wall: uniform source ~ 0 at this setup.
    let mut turb = turb_state(&g, 3.0, 1e6, 0.01);
    advance_turb(&mut turb, &st, &g, &bc, 1e-4).unwrap();
    for x in &turb.nu_tilde {
        assert!((x - 3.0).abs() < 1e-9, "got {x}");
    }
}

#[test]
fn zero_nu_tilde_is_a_fixed_point() {
    // the laminar state nu_tilde = 0 must stay exactly zero: production is
    // linear in nu_tilde, destruction quadratic, so zero is preserved.
    let g = grid2d(10, 10, 0.0, 1.0, 0.0, 1.0).unwrap();
    let st = uniform_state(&g, 1.0, 0.3, 0.2, 1.0);
    let bc = Boundaries2d::default();
    let mut turb = turb_state(&g, 0.0, 0.05, 0.01);
    advance_turb(&mut turb, &st, &g, &bc, 1e-4).unwrap();
    for x in &turb.nu_tilde {
        assert!(x.abs() < 1e-10, "got {x}");
    }
}

#[test]
fn wall_proximity_drives_nu_tilde_down() {
    // near a wall (small d) with zero vorticity, destruction dominates
    // production, so nu_tilde decreases.
    let g = grid2d(8, 8, 0.0, 1.0, 0.0, 1.0).unwrap();
    let st = uniform_state(&g, 1.0, 0.3, 0.0, 1.0);
    let bc = Boundaries2d {
        south: Bc2d::NoSlipWall,
        ..Default::default()
    };
    let n = g.nx * g.ny;
    let mut turb = TurbState {
        nu_tilde: vec![3.0; n],
        d: g.centers_y.iter().map(|y| (*y).max(1e-3)).collect(),
        params: SaParams {
            mu: 0.01,
            pr: 0.72,
            nu_tilde_inf: 3.0,
            pr_t: 0.9,
        },
    };
    advance_turb(&mut turb, &st, &g, &bc, 1e-3).unwrap();
    let near = turb.nu_tilde[0]; // bottom row, closest to the "wall"
    let far = turb.nu_tilde[g.nx * (g.ny - 1)];
    assert!(near < 3.0, "near-wall nu_tilde must drop, got {near}");
    assert!(far <= 3.0 + 1e-9, "far cells must not grow, got {far}");
}

#[test]
fn advection_carries_a_front_left_to_right() {
    // a left half at 5, right half at 1, uniform u > 0: after one step the
    // front must move right (the upwind scheme transports nu_tilde downstream).
    let g = grid2d(16, 4, 0.0, 1.0, 0.0, 0.25).unwrap();
    let st = uniform_state(&g, 1.0, 0.5, 0.0, 1.0);
    let n = g.nx * g.ny;
    let mut nt = vec![1.0; n];
    for j in 0..g.ny {
        for i in 0..g.nx / 2 {
            nt[j * g.nx + i] = 5.0;
        }
    }
    let bc = Boundaries2d::default();
    let mut turb = TurbState {
        nu_tilde: nt,
        d: vec![1e6; n],
        params: SaParams {
            mu: 0.01,
            pr: 0.72,
            nu_tilde_inf: 1.0,
            pr_t: 0.9,
        },
    };
    advance_turb(&mut turb, &st, &g, &bc, 1e-3).unwrap();
    let front = g.nx / 2;
    let mid = g.ny / 2;
    let after_front = turb.nu_tilde[mid * g.nx + front];
    let before_front = turb.nu_tilde[mid * g.nx + front - 1];
    assert!(
        after_front > 1.0,
        "front cell must gain nu_tilde, got {after_front}"
    );
    assert!(
        before_front < 5.0,
        "behind-front cell must lose, got {before_front}"
    );
}

#[test]
fn wall_distance_is_perpendicular_to_solid_sides() {
    // south wall only: d equals the cell-center y exactly; other cells grow
    // with y; a no-wall domain gets the far value everywhere.
    let g = grid2d(4, 8, 0.0, 1.0, 0.0, 1.0).unwrap();
    let bc_wall = Boundaries2d {
        south: Bc2d::NoSlipWall,
        ..Default::default()
    };
    let d = wall_distance2d(&g, &bc_wall, 100.0);
    for j in 0..g.ny {
        for i in 0..g.nx {
            let c = j * g.nx + i;
            assert!((d[c] - g.centers_y[c]).abs() < 1e-12);
        }
    }
    let bc_open = Boundaries2d::default();
    let d_open = wall_distance2d(&g, &bc_open, 100.0);
    assert!(d_open.iter().all(|x| (*x - 100.0).abs() < 1e-12));
}

#[test]
fn diffusion_decays_a_peak_at_the_analytic_rate() {
    // a gaussian peak in nu_tilde with zero velocity: the active term is the
    // (nu + nu_tilde)/sigma laplacian. the reference is the gaussian moment
    // closure (exact for constant diffusivity, a good closure when the sa
    // term makes diffusivity amplitude-dependent): dA/dt = -2 D A / s2, with
    // d s2/dt = 2 D and D = (nu + A)/sigma. dt respects the explicit bound
    // dt < dx^2 / (2 D_max).
    use nufor_core::SIGMA;
    let nx = 64;
    let g = grid2d(nx, 4, 0.0, 1.0, 0.0, 0.0625).unwrap();
    let st = uniform_state(&g, 1.0, 0.0, 0.0, 1.0);
    let mid = g.ny / 2;
    let x0 = 0.5;
    let sig0 = 0.06;
    let amp = 0.05;
    let mut nt = vec![0.0f64; g.nx * g.ny];
    for i in 0..g.nx {
        let x = g.centers_x[mid * g.nx + i];
        let val = amp * (-(x - x0) * (x - x0) / (2.0 * sig0 * sig0)).exp();
        for j in 0..g.ny {
            nt[j * g.nx + i] = val;
        }
    }
    let bc = Boundaries2d::default();
    let nu = 0.5;
    let mut turb = TurbState {
        nu_tilde: nt,
        d: vec![1e6; g.nx * g.ny],
        params: SaParams {
            mu: nu,
            pr: 0.72,
            nu_tilde_inf: amp,
            pr_t: 0.9,
        },
    };
    let dt = 3e-5;
    let steps = 200;
    for _ in 0..steps {
        advance_turb(&mut turb, &st, &g, &bc, dt).unwrap();
    }
    // constant-diffusivity reference (nu_tilde << nu so d ~ nu/sigma):
    // sigma^2(t) = sigma0^2 + 2 d t, peak = amp * sigma0 / sigma(t).
    let dcoef = nu / SIGMA;
    let t = steps as f64 * dt;
    let sig = (sig0 * sig0 + 2.0 * dcoef * t).sqrt();
    let peak_ref = amp * sig0 / sig;
    let peak_num = turb.nu_tilde[mid * g.nx + nx / 2];
    let rel = (peak_num - peak_ref).abs() / peak_ref;
    assert!(rel < 0.05, "peak {peak_num} vs ref {peak_ref}, rel {rel}");
}

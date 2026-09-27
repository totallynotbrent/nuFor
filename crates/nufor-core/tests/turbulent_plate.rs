//! the turbulent flat-plate acceptance: the sa model carrying a genuinely
//! turbulent layer at correlation reynolds numbers. a short march on a
//! coarse-but-honest grid asserts the layer develops (nu_tilde grows well
//! past molecular), the wall shear tracks the power-law correlation within
//! the band the full-resolution march measured, and nothing degenerates.
//!
//! the march is a real boundary-layer integration (tens of thousands of
//! explicit steps even at this reduced size), far too slow for the ci
//! gate's debug-mode suite, so the test is marked ignore and runs
//! explicitly with `cargo test --release -- --ignored`.

use nufor_core::{
    advance2d_sa_rk2, cons_to_prim2d, prim_to_cons2d, stretched_grid2d, wall_distance2d, Bc2d,
    BlasiusProfile, Boundaries2d, Clustering, ConservedState2d, InflowProfile, SaParams, TurbState,
};

const GAMMA: f64 = 1.4;
const U_INF: f64 = 0.2;
const MU: f64 = 4.8e-7;
const RHO: f64 = 1.0;
const P_INF: f64 = 1.0;

/// the discrete cf at every station: the quadratic-consistent wall
/// derivative the solver applies, reported against the freestream head.
fn wall_cf(g: &nufor_core::Grid2d, st: &ConservedState2d) -> Vec<f64> {
    let (u, _, _) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    let mut cf = vec![0.0; g.nx];
    for i in 0..g.nx {
        let (y0, y1) = (g.centers_y[i], g.centers_y[g.nx + i]);
        let y_wall = g.ymin;
        let (u0, u1) = (u[i], u[g.nx + i]);
        let w0 = (y_wall - y1) / ((y0 - y_wall) * (y0 - y1));
        let w1 = (y_wall - y0) / ((y1 - y_wall) * (y1 - y0));
        let du = u0 * w0 + u1 * w1;
        cf[i] = 2.0 * MU * du / (RHO * U_INF * U_INF);
    }
    cf
}

#[test]
#[ignore = "a full boundary-layer march; run with --release -- --ignored"]
fn sa_grows_a_turbulent_layer_with_correlation_cf() {
    // the full-resolution march (160x60, t=4) measured: nu_tilde reaching
    // hundreds of nu at mid-chord and mid-domain cf 12-19 percent above
    // the schlichting power law. this test runs the same physics at a
    // quarter resolution and a shorter march, asserting the same
    // contracts with margins widened to cover the coarser grid: the layer
    // is turbulent (nu_tilde >> nu), the cf level is in the correlation's
    // neighborhood (within 60 percent), and the state stays clean.
    let nx = 40;
    let ny = 30;
    let g = stretched_grid2d(
        nx,
        0.0,
        1.2,
        ny,
        0.0,
        0.35,
        Clustering {
            first_cell: 5.18e-5,
            growth: 1.2,
        },
    )
    .unwrap();
    let n = nx * ny;
    let nu = MU / RHO;
    let nu_tilde_inf = 3.0 * nu;
    let et = vec![P_INF / (GAMMA - 1.0) + 0.5 * U_INF * U_INF; n];
    let (mx, my, e) = prim_to_cons2d(&vec![RHO; n], &vec![U_INF; n], &vec![0.0; n], &et).unwrap();
    let mut st = ConservedState2d {
        rho: vec![RHO; n],
        mx,
        my,
        e,
    };
    // the leading edge sits 0.2 upstream: the layer is singular at the
    // edge itself, and the inflow plane must be downstream of it.
    let p = BlasiusProfile::new(U_INF, nu, -0.2)
        .unwrap()
        .anchored_at(0.0)
        .unwrap();
    let bc = Boundaries2d {
        west: Bc2d::ProfileInflow {
            profile: std::sync::Arc::new(InflowProfile::Blasius(p)),
            rho: RHO,
            p: P_INF,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::NoSlipWall,
        north: Bc2d::SupersonicOutflow,
    };
    let mut turb = TurbState {
        nu_tilde: vec![nu_tilde_inf; n],
        d: wall_distance2d(&g, &bc, 0.35),
        params: SaParams {
            mu: MU,
            pr: 0.72,
            pr_t: 0.9,
            nu_tilde_inf,
        },
    };

    // march one flow-through (t = 1.2 / u_inf = 6): enough for the layer
    // to develop past transition at this reynolds range.
    let mut t = 0.0f64;
    while t < 6.0 {
        t += advance2d_sa_rk2(&mut st, &mut turb, &g, GAMMA, 0.3, true, &bc).unwrap();
    }

    // the layer must be turbulent: nu_tilde grows far past molecular
    // viscosity at mid-chord (the march measured ~220 nu at full grid).
    let mid = n / 2;
    let nt_mid = turb.nu_tilde[mid];
    assert!(
        nt_mid > 10.0 * nu,
        "nu_tilde at mid-chord {nt_mid:.3e} must exceed 10 nu ({:.3e})",
        10.0 * nu
    );

    // the wall shear must be in the correlation's neighborhood over the
    // interior band. the corners carry the bounded-domain compatibility
    // transient (the laminar hold documents the same artifact), and on
    // this coarse grid it reaches further in: only the middle half of the
    // plate is a clean comparison.
    let cf = wall_cf(&g, &st);
    let band = nx / 4..3 * nx / 4;
    let mut worst = 0.0f64;
    for i in band {
        let re_x = RHO * U_INF * g.centers_x[i] / MU;
        let cf_corr = 0.0592 * re_x.powf(-0.2);
        let ratio = cf[i] / cf_corr;
        worst = worst.max((ratio - 1.0).abs());
        assert!(cf[i] > 0.0, "cf must be positive at station {i}");
    }
    assert!(
        worst < 0.6,
        "worst correlation deviation {worst:.3} beyond the 60 percent band"
    );

    // and the state must be clean.
    assert!(st.rho.iter().all(|r| r.is_finite() && *r > 0.0));
    assert!(turb.nu_tilde.iter().all(|x| x.is_finite() && *x >= 0.0));
}

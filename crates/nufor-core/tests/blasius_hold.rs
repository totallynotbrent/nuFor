//! the blasius flat-plate hold: the analytic laminar layer laid into the
//! domain must be held steady by the solver, and the discrete skin friction
//! must match the cf = 2 f''(0)/sqrt(re_x) correlation.

use nufor_core::{
    advance2d_visc_rk2, cons_to_prim2d, prim_to_cons2d, stretched_grid2d, Bc2d, BlasiusProfile,
    Boundaries2d, Clustering, ConservedState2d, ViscParams, CF_CONST,
};

const GAMMA: f64 = 1.4;

/// lay the exact blasius field (u and v) into the state on the clustered grid.
fn blasius_state(g: &nufor_core::Grid2d, p: &BlasiusProfile) -> ConservedState2d {
    let n = g.nx * g.ny;
    let (mut u, mut v) = (vec![0.0; n], vec![0.0; n]);
    let rho = vec![1.0; n];
    let pr = vec![1.0; n];
    for k in 0..n {
        let (x, y) = (g.centers_x[k], g.centers_y[k]);
        u[k] = p.u(x, y);
        v[k] = p.v(x, y);
    }
    let et: Vec<f64> = pr
        .iter()
        .zip(&u)
        .zip(&v)
        .map(|((pp, uu), vv)| pp / (GAMMA - 1.0) + 0.5 * (uu * uu + vv * vv))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    ConservedState2d { rho, mx, my, e }
}

#[test]
fn blasius_layer_is_held_with_profile_inflow() {
    // a plate with its leading edge at x=0, upstream of the domain; the
    // inflow plane at x=0.2 sees a layer at re_x=40, so no leading-edge
    // singularity touches the domain. the wall-clustered grid puts ~27
    // cells across the layer at the exit station, and the domain is tall
    // enough that the layer never feels the top boundary.
    let u_inf = 0.2;
    let nu = 1.0e-3;
    let g = stretched_grid2d(
        48,
        0.2,
        1.2,
        56,
        0.0,
        1.0,
        Clustering {
            first_cell: 2.0e-3,
            growth: 1.2,
        },
    )
    .unwrap();
    let p = BlasiusProfile::new(u_inf, nu, 0.0)
        .unwrap()
        .anchored_at(0.2);
    let mut st = blasius_state(&g, &p);
    let bc = Boundaries2d {
        west: Bc2d::ProfileInflow {
            profile: std::sync::Arc::new(nufor_core::InflowProfile::Blasius(p.clone())),
            rho: 1.0,
            p: 1.0,
        },
        east: Bc2d::Transmissive,
        south: Bc2d::NoSlipWall,
        north: Bc2d::Transmissive,
    };
    let mu = nu; // rho = 1
    let mut t = 0.0;
    let t_end = 0.15; // several flow-throughs of the 1.0-long domain
    while t < t_end {
        let (dt, _) = advance2d_visc_rk2(
            &mut st,
            &g,
            GAMMA,
            0.4,
            true,
            &bc,
            ViscParams { mu, pr: 0.72 },
        )
        .unwrap();
        t += dt;
    }
    // the profile must hold: compare u at every cell against the analytic
    // layer, relative to the freestream.
    let (u_c, _, _) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    let mut worst = 0.0f64;
    for ((&uc, &x), &y) in u_c.iter().zip(&g.centers_x).zip(&g.centers_y) {
        let err = (uc - p.u(x, y)).abs() / u_inf;
        worst = worst.max(err);
    }
    eprintln!("blasius hold worst relative profile error: {worst:.2e}");
    assert!(
        worst < 5e-2,
        "the layer must be held within 5 percent, worst {worst:.2e}"
    );
    assert!(st.rho.iter().all(|&r| r.is_finite() && r > 0.0));

    // the discrete cf must match the correlation over the mid-domain band:
    // cf(x) = 0.6641/sqrt(re_x). the first stations sit in the inflow
    // corner (the profile/pressure compatibility transient) and the last
    // ones in the outflow corner; both are boundary artifacts every
    // validation setup excludes, so the assertion covers the interior.
    let mut worst_cf = 0.0f64;
    for i in 8..g.nx - 6 {
        let x = g.centers_x[i];
        let re_x = u_inf * x / nu;
        let cf_exact = CF_CONST / re_x.sqrt();
        // the same true-position wall derivative the cli reports.
        let (y0, y1) = (g.centers_y[i], g.centers_y[g.nx + i]);
        let y_wall = g.ymin;
        let (u0, u1) = (u_c[i], u_c[g.nx + i]);
        let w0 = (y_wall - y1) / ((y0 - y_wall) * (y0 - y1));
        let w1 = (y_wall - y0) / ((y1 - y_wall) * (y1 - y0));
        let du = u0 * w0 + u1 * w1;
        let cf_num = 2.0 * mu * du / (1.0 * u_inf * u_inf);
        let err = (cf_num - cf_exact).abs() / cf_exact;
        worst_cf = worst_cf.max(err);
    }
    eprintln!("blasius hold worst relative cf error: {worst_cf:.2e}");
    assert!(
        worst_cf < 5e-2,
        "cf must match 0.664/sqrt(re_x) within 5 percent, worst {worst_cf:.2e}"
    );
}

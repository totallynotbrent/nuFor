//! stretched-grid solver tests: the exact analytic holds must survive
//! wall-normal clustering, the uniform fast path must stay consistent with
//! the plain grid, and the generators must produce sane face coordinates.

use nufor_core::{
    advance2d, advance2d_visc_rk2, channel_grid2d, cons_to_prim2d, prim_to_cons2d,
    rectilinear_grid2d, stretched_grid2d, Bc2d, Boundaries2d, Clustering, ConservedState2d,
    ViscParams,
};

const GAMMA: f64 = 1.4;

#[test]
fn stretched_grid_face_coordinates_are_sane() {
    // growth 1 with h0 = uniform height reproduces plain spacing.
    let g = stretched_grid2d(
        8,
        0.0,
        1.0,
        8,
        0.0,
        1.0,
        Clustering {
            first_cell: 0.125,
            growth: 1.0,
        },
    )
    .unwrap();
    assert!((g.faces_x[1] - g.faces_x[0]) - 0.125 < 1e-12);
    assert!((g.faces_y[1] - g.faces_y[0]) - 0.125 < 1e-12);
    // growth > 1 clusters near the lower wall: first cell the smallest.
    let s = stretched_grid2d(
        8,
        0.0,
        1.0,
        8,
        0.0,
        1.0,
        Clustering {
            first_cell: 0.01,
            growth: 1.2,
        },
    )
    .unwrap();
    assert!(s.dys[0] < s.dys[1]);
    assert!(s.dys.iter().all(|w| *w > 0.0));
    assert!((s.faces_y.last().unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn channel_grid_clusters_both_walls() {
    let g = channel_grid2d(
        12,
        0.0,
        0.5,
        40,
        0.0,
        1.0,
        Clustering {
            first_cell: 0.01,
            growth: 1.15,
        },
    )
    .unwrap();
    assert!(g.dys[0] < g.dys[g.ny / 2]);
    assert!(g.dys[g.ny - 1] < g.dys[g.ny / 2]);
    assert!((g.faces_y.last().unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn poiseuille_parabola_is_held_on_a_stretched_channel() {
    // channel [0, 0.5] x [0, 1] with both walls clustered; the analytic
    // parabola must still be held steady against the body force.
    let g = channel_grid2d(
        24,
        0.0,
        0.5,
        40,
        0.0,
        1.0,
        Clustering {
            first_cell: 0.008,
            growth: 1.15,
        },
    )
    .unwrap();
    let mu = 0.05;
    let umax = 0.05;
    let body = 8.0 * mu * umax;
    let n = g.nx * g.ny;
    let mut u = vec![0.0; n];
    for j in 0..g.ny {
        let y = g.centers_y[j * g.nx];
        for i in 0..g.nx {
            u[j * g.nx + i] = umax * 4.0 * y * (1.0 - y);
        }
    }
    let rho = vec![1.0; n];
    let v = vec![0.0; n];
    let p = vec![1.0; n];
    let et: Vec<f64> = p
        .iter()
        .zip(&u)
        .map(|(pp, uu)| *pp / (GAMMA - 1.0) + 0.5 * uu * uu)
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let bc = Boundaries2d {
        south: Bc2d::NoSlipWall,
        north: Bc2d::NoSlipWall,
        west: Bc2d::Transmissive,
        east: Bc2d::Transmissive,
    };
    let mut t = 0.0;
    while t < 0.01 {
        let (dt, _) = advance2d_visc_rk2(
            &mut st,
            &g,
            GAMMA,
            0.5,
            true,
            &bc,
            ViscParams { mu, pr: 0.72 },
        )
        .unwrap();
        for mxv in st.mx.iter_mut() {
            *mxv += dt * body;
        }
        t += dt;
    }
    let (u_c, _, _) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    let mid = g.nx / 2;
    let mut worst = 0.0f64;
    for j in 0..g.ny {
        let y = g.centers_y[j * g.nx];
        let err = (u_c[j * g.nx + mid] - umax * 4.0 * y * (1.0 - y)).abs();
        worst = worst.max(err);
    }
    eprintln!("stretched poiseuille worst error: {worst:.2e}");
    assert!(
        worst < 1e-2,
        "stretched channel must hold the parabola, worst {worst:.2e}"
    );
    assert!(st.rho.iter().all(|&r| r.is_finite() && r > 0.0));
}

#[test]
fn uniform_grid_is_the_unchanged_fast_path() {
    // a uniform grid built through the stretched generator must equal the
    // plain grid2d on every metric the solvers touch.
    let s = stretched_grid2d(
        16,
        0.0,
        2.0,
        16,
        0.0,
        2.0,
        Clustering {
            first_cell: 0.125,
            growth: 1.0,
        },
    )
    .unwrap();
    let p = nufor_core::grid2d(16, 16, 0.0, 2.0, 0.0, 2.0).unwrap();
    assert_eq!(s.nx, p.nx);
    assert_eq!(s.ny, p.ny);
    for k in 0..s.nx * s.ny {
        assert!((s.centers_x[k] - p.centers_x[k]).abs() < 1e-12);
        assert!((s.centers_y[k] - p.centers_y[k]).abs() < 1e-12);
    }
    for i in 0..=s.nx {
        assert!((s.faces_x[i] - p.faces_x[i]).abs() < 1e-12);
    }
    for j in 0..=s.ny {
        assert!((s.faces_y[j] - p.faces_y[j]).abs() < 1e-12);
    }
}

#[test]
fn rectilinear_grid2d_accepts_stretched_axes() {
    // the mesh-import path: per-axis face coordinates with real stretching.
    let mut fy = vec![0.0f64];
    let mut y = 0.0f64;
    let mut w = 0.01f64;
    while y < 1.0 {
        y += w;
        fy.push(y.min(1.0));
        w *= 1.2;
    }
    let fx: Vec<f64> = (0..=10).map(|i| i as f64 * 0.1).collect();
    let g = rectilinear_grid2d(&fx, &fy).unwrap();
    assert!(g.dys[0] < g.dys[g.dys.len() - 1]);
    assert!(g.dys.iter().all(|w| *w > 0.0));
}

#[test]
fn euler_step_is_stable_and_conservative_on_a_stretched_grid() {
    // a gaussian bump advected on a one-sided clustered grid, small enough
    // that no wave reaches the boundary: the state stays finite and the
    // volume-weighted mass holds to roundoff (transmissive sides telescope).
    let g = stretched_grid2d(
        60,
        0.0,
        3.0,
        60,
        0.0,
        3.0,
        Clustering {
            first_cell: 0.02,
            growth: 1.15,
        },
    )
    .unwrap();
    let n = g.nx * g.ny;
    let rho: Vec<f64> = (0..n)
        .map(|k| {
            let (x, y) = (g.centers_x[k], g.centers_y[k]);
            1.0 + 0.5 * (-8.0 * ((x - 1.5).powi(2) + (y - 1.5).powi(2))).exp()
        })
        .collect();
    let u = vec![0.5; n];
    let v = vec![0.3; n];
    let p = vec![1.0; n];
    let et: Vec<f64> = p
        .iter()
        .zip(&u)
        .zip(&v)
        .map(|((pp, uu), vv)| *pp / (GAMMA - 1.0) + 0.5 * (uu * uu + vv * vv))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let bc = Boundaries2d::default();
    let vol: Vec<f64> = (0..n)
        .map(|k| {
            let (i, j) = (k % g.nx, k / g.nx);
            g.dxs[i] * g.dys[j]
        })
        .collect();
    let m0: f64 = st.rho.iter().zip(&vol).map(|(r, w)| r * w).sum();
    for _ in 0..40 {
        advance2d(&mut st, &g, GAMMA, 0.4, true, &bc).unwrap();
    }
    let m1: f64 = st.rho.iter().zip(&vol).map(|(r, w)| r * w).sum();
    assert!(st.rho.iter().all(|&r| r.is_finite() && r > 0.0));
    let drift = (m1 - m0).abs() / m0;
    assert!(drift < 1e-6, "mass drifted {drift} on a stretched grid");
}

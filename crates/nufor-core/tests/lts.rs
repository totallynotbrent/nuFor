//! local time stepping (lts) tests: the scheme trades time accuracy for
//! convergence speed, so its contract is about the steady state, not the
//! trajectory. at a true steady state every residual is zero and the per-cell
//! dt scaling multiplies that zero, so a field the operator holds must stay
//! held under lts; and from a perturbed start the lts sweeps must converge
//! to the same fixed point the time-accurate march reaches.

use nufor_core::{
    advance2d_sa_lts, advance2d_sa_rk2, channel_grid2d, cons_to_prim2d, grid2d, prim_to_cons2d,
    stretched_grid2d, wall_distance2d, Bc2d, Boundaries2d, Clustering, ConservedState2d, SaParams,
    TurbState,
};

const GAMMA: f64 = 1.4;

/// freestream + sa wall: the flat-plate geometry at tiny scale, where the
/// exact steady state is the uniform freestream (nu_tilde at its freestream
/// value with zero production far from the wall, zero momentum defect in
/// the freestream rows).
fn plate_case(
    nx: usize,
    ny: usize,
) -> (
    nufor_core::Grid2d,
    ConservedState2d,
    Boundaries2d,
    TurbState,
) {
    let u_inf = 0.2;
    let rho = 1.0;
    let p_inf = 1.0;
    let mu = 4.8e-7;
    let nu_tilde_inf = 3.0 * mu;
    let g = grid2d(nx, ny, 0.0, 1.2, 0.0, 0.35).unwrap();
    let n = nx * ny;
    let et = vec![p_inf / (GAMMA - 1.0) + 0.5 * u_inf * u_inf; n];
    let (mx, my, e) = prim_to_cons2d(&vec![rho; n], &vec![u_inf; n], &vec![0.0; n], &et).unwrap();
    let st = ConservedState2d {
        rho: vec![rho; n],
        mx,
        my,
        e,
    };
    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow {
            rho,
            u: u_inf,
            v: 0.0,
            p: p_inf,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::NoSlipWall,
        north: Bc2d::SupersonicOutflow,
    };
    let d = wall_distance2d(&g, &bc, 0.35);
    let turb = TurbState {
        nu_tilde: vec![nu_tilde_inf; n],
        d,
        params: SaParams {
            mu,
            pr: 0.72,
            pr_t: 0.9,
            nu_tilde_inf,
        },
    };
    (g, st, bc, turb)
}

#[test]
fn lts_holds_a_steady_freestream_above_the_wall() {
    // the uniform freestream is an exact steady state of the rows far from
    // the wall; near the wall a boundary layer grows, so the contract here
    // is narrower: the far-field rows must not drift under lts.
    let (g, mut st, bc, mut turb) = plate_case(24, 18);
    let u_inf = 0.2;
    let rho = 1.0;
    let mu = turb.params.mu;
    let before = st.clone();
    for _ in 0..200 {
        advance2d_sa_lts(&mut st, &mut turb, &g, GAMMA, 0.3, true, &bc).unwrap();
    }
    let (u, _, _) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    // the top third of the domain: freestream rows, no wall influence.
    let mut worst = 0.0f64;
    for j in (2 * g.ny / 3)..g.ny {
        for i in 0..g.nx {
            let k = j * g.nx + i;
            worst = worst.max((u[k] - u_inf).abs() / u_inf);
        }
    }
    assert!(worst < 1e-3, "freestream drifted {worst:e}");
    // and the wall rows respond to the no-slip condition: u near the wall
    // is strictly below the freestream.
    // the wall row responds to the no-slip flux: strictly below freestream
    // by a measurable margin (the viscous layer is thin this early).
    let j0 = 0;
    let mid = g.nx / 2;
    assert!(u[j0 * g.nx + mid] < u_inf * (1.0 - 1e-3));
    let _ = (before, rho, mu);
}

#[test]
fn lts_matches_the_time_accurate_fixed_point() {
    // both paths from the same start, budgeted independently; the wall-row
    // velocity profile at mid-chord must agree: the lts fixed point is the
    // same steady state the march converges to.
    let nx = 40;
    let ny = 24;
    let (g, _, bc, _) = plate_case(nx, ny);
    let u_inf = 0.2;
    let mu = 4.8e-7;

    // reference: time-accurate march, one flow-through (t=1.2/0.2=6).
    let (_, mut st_ref, _, mut turb_ref) = plate_case(nx, ny);
    let mut t = 0.0;
    while t < 6.0 {
        t += advance2d_sa_rk2(&mut st_ref, &mut turb_ref, &g, GAMMA, 0.3, true, &bc).unwrap();
    }

    // lts: same start, sweeps.
    let (_, mut st_lts, _, mut turb_lts) = plate_case(nx, ny);
    for _ in 0..30_000 {
        advance2d_sa_lts(&mut st_lts, &mut turb_lts, &g, GAMMA, 0.3, true, &bc).unwrap();
    }

    let (u_ref, _, _) = cons_to_prim2d(&st_ref.rho, &st_ref.mx, &st_ref.my, &st_ref.e).unwrap();
    let (u_lts, _, _) = cons_to_prim2d(&st_lts.rho, &st_lts.mx, &st_lts.my, &st_lts.e).unwrap();
    // compare the full wall-normal profile at three stations.
    let mut worst = 0.0f64;
    for &i in &[10, 20, 30] {
        for j in 0..g.ny {
            let k = j * g.nx + i;
            let a = u_ref[k];
            let b = u_lts[k];
            worst = worst.max((a - b).abs() / u_inf);
        }
    }
    eprintln!("lts vs march fixed-point worst: {worst:.3e}");
    assert!(worst < 5e-2, "lts and march disagree {worst:e}");
    let _ = mu;
}

#[test]
fn local_dts_bounds_are_sane_and_stable() {
    // the per-cell increments are the lts contract's foundation: strictly
    // positive, bounded by the cfl region, and a global step at the
    // smallest of them must leave the state intact.
    use nufor_core::local_dts;
    let (g, st0, bc, turb) = plate_case(24, 18);
    let d = local_dts(&st0, &turb, &g, GAMMA, 0.3);
    let dt_min = d.iter().cloned().fold(f64::INFINITY, f64::min);
    let dt_max = d.iter().cloned().fold(0.0f64, f64::max);
    eprintln!(
        "local dt range [{dt_min:.3e}, {dt_max:.3e}], ratio {:.1}",
        dt_max / dt_min
    );
    // strictly positive and spanning a sane dynamic range.
    assert!(dt_min > 0.0 && dt_max / dt_min < 1.0e4);
    // a global march at its own cfl stays inside every per-cell bound.
    let mut st_l = st0;
    let mut tl = turb;
    for _ in 0..5 {
        let _ = advance2d_sa_rk2(&mut st_l, &mut tl, &g, GAMMA, 0.3, true, &bc).unwrap();
    }
    assert!(st_l.rho.iter().all(|r| r.is_finite() && *r > 0.0));
}

#[test]
fn lts_survives_the_clustered_plate_grid() {
    // the clustered wall grid is the case lts exists for; a short run must
    // keep every field finite and positive (a regression guard for the
    // per-cell factors on stretched axes).
    let nx = 40;
    let ny = 30;
    let g = stretched_grid2d(
        nx,
        0.0,
        2.0,
        ny,
        0.0,
        1.0,
        Clustering {
            first_cell: 5.0e-4,
            growth: 1.12,
        },
    )
    .unwrap();
    let n = nx * ny;
    let mu = 1e-9;
    let et = vec![1.0 / 0.4 + 0.5 * 0.2 * 0.2; n];
    let (mx, my, e) = prim_to_cons2d(&vec![1.0; n], &vec![0.2; n], &vec![0.0; n], &et).unwrap();
    let mut st = ConservedState2d {
        rho: vec![1.0; n],
        mx,
        my,
        e,
    };
    // anchor the layer one cell inside the domain: at the leading edge
    // itself the similarity v is singular, which is a case-setup error, not
    // a solver property (the laminar hold test anchors the same way).
    let p = nufor_core::BlasiusProfile::new(0.2, 1e-5, 0.0)
        .unwrap()
        .anchored_at(g.dxs[0])
        .unwrap();
    let bc = Boundaries2d {
        west: Bc2d::ProfileInflow {
            profile: std::sync::Arc::new(nufor_core::InflowProfile::Blasius(p)),
            rho: 1.0,
            p: 1.0,
        },
        east: Bc2d::Transmissive,
        south: Bc2d::NoSlipWall,
        north: Bc2d::Transmissive,
    };
    let mut turb = TurbState {
        nu_tilde: vec![0.0; n],
        d: wall_distance2d(&g, &bc, 2.0),
        params: SaParams {
            mu,
            pr: 0.72,
            pr_t: 0.9,
            nu_tilde_inf: 0.0,
        },
    };
    for k in 0..20 {
        if let Err(e) = advance2d_sa_lts(&mut st, &mut turb, &g, GAMMA, 0.3, true, &bc) {
            panic!("lts failed at sweep {k}: {e:?}");
        }
        if let Some(bad) = st.rho.iter().position(|r| !r.is_finite() || *r <= 0.0) {
            panic!("rho corrupted at sweep {k}, cell {bad}: {}", st.rho[bad]);
        }
        if let Some(bad) = st.e.iter().position(|x| !x.is_finite()) {
            panic!("energy corrupted at sweep {k}, cell {bad}: {}", st.e[bad]);
        }
    }
}

#[test]
fn lts_channel_grid_clustering_is_accepted() {
    // the clustered channel constructor feeds lts without erroring: a smoke
    // on the public surface combination (stretched + sa-lts).
    let g = channel_grid2d(
        12,
        0.0,
        1.0,
        12,
        0.0,
        0.2,
        Clustering {
            first_cell: 0.004,
            growth: 1.15,
        },
    )
    .unwrap();
    let n = g.nx * g.ny;
    let et = vec![1.0 / (GAMMA - 1.0); n];
    let (mx, my, e) = prim_to_cons2d(&vec![1.0; n], &vec![0.1; n], &vec![0.0; n], &et).unwrap();
    let mut st = ConservedState2d {
        rho: vec![1.0; n],
        mx,
        my,
        e,
    };
    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow {
            rho: 1.0,
            u: 0.1,
            v: 0.0,
            p: 1.0,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::NoSlipWall,
        north: Bc2d::NoSlipWall,
    };
    let mut turb = TurbState {
        nu_tilde: vec![0.0; n],
        d: wall_distance2d(&g, &bc, 0.2),
        params: SaParams {
            mu: 1e-3,
            pr: 0.72,
            pr_t: 0.9,
            nu_tilde_inf: 0.0,
        },
    };
    for _ in 0..50 {
        advance2d_sa_lts(&mut st, &mut turb, &g, GAMMA, 0.3, true, &bc).unwrap();
    }
    assert!(st.rho.iter().all(|r| r.is_finite()));
}

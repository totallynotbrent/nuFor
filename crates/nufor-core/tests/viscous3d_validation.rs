//! the 3d validation battery for the viscous operator: poiseuille channel
//! flow held steady, and the wall-brake behavior end-to-end through the
//! heun march.

use nufor_core::{
    advance3d_visc_rk2, cons_to_prim3d, grid3d, prim_to_cons3d, Bounds3d, ConservedState3d, Grid3d,
    ViscParams3d, Walls3d,
};

const GAMMA: f64 = 1.4;

/// the parabolic poiseuille field u(y) between plates at y=0 and y=1.
fn parabola(g: &Grid3d, umax: f64) -> Vec<f64> {
    let n = g.nx * g.ny * g.nz;
    let mut u = vec![0.0; n];
    for k in 0..g.nz {
        for j in 0..g.ny {
            let y = g.centers_y[(k * g.ny + j) * g.nx];
            for i in 0..g.nx {
                u[(k * g.ny + j) * g.nx + i] = umax * 4.0 * y * (1.0 - y);
            }
        }
    }
    u
}

#[test]
fn poiseuille_3d_parabola_is_held_steady() {
    // a slab of channel [0, 0.5] x [0, 1] x [0, 0.25]: no-slip plates at
    // y=0 and y=1, transmissive everywhere else (spanwise-invariant flow).
    // the body force B = 8 mu umax / H^2 balances mu d2u/dy2 so the
    // parabola is steady; the march must hold it.
    let b = Bounds3d {
        xmin: 0.0,
        xmax: 0.5,
        ymin: 0.0,
        ymax: 1.0,
        zmin: 0.0,
        zmax: 0.25,
    };
    let g = grid3d(16, 24, 6, &b).unwrap();
    let (mu, umax) = (0.05, 0.05);
    let body = 8.0 * mu * umax;
    let u = parabola(&g, umax);
    let n = g.nx * g.ny * g.nz;
    let rho = vec![1.0; n];
    let v = vec![0.0; n];
    let w = vec![0.0; n];
    let et: Vec<f64> = u
        .iter()
        .map(|uu| 1.0 / (GAMMA - 1.0) + 0.5 * uu * uu)
        .collect();
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &u, &v, &w, &et).unwrap();
    let mut st = ConservedState3d { rho, mx, my, mz, e };
    let walls = Walls3d {
        south: true,
        north: true,
        ..Default::default()
    };
    let mut t = 0.0;
    while t < 0.05 {
        let (dt, _) = advance3d_visc_rk2(
            &mut st,
            &g,
            GAMMA,
            0.5,
            true,
            ViscParams3d { mu, pr: 0.72 },
            &walls,
        )
        .unwrap();
        for mxv in st.mx.iter_mut() {
            *mxv += dt * body;
        }
        t += dt;
    }
    // the profile must stay the parabola at the channel mid-plane and at
    // every spanwise station (spanwise invariance).
    let (u_c, _, _, _) = cons_to_prim3d(&st.rho, &st.mx, &st.my, &st.mz, &st.e).unwrap();
    let mid = g.nx / 2;
    let mut worst = 0.0f64;
    for k in 0..g.nz {
        for j in 0..g.ny {
            let y = g.centers_y[(k * g.ny + j) * g.nx + mid];
            let err = (u_c[(k * g.ny + j) * g.nx + mid] - umax * 4.0 * y * (1.0 - y)).abs();
            worst = worst.max(err);
        }
    }
    eprintln!("3d poiseuille worst profile error: {worst:.2e}");
    assert!(
        worst < 1e-2,
        "3d channel must hold the parabola, worst error {worst:.2e}"
    );
    assert!(
        st.rho.iter().all(|&r| r.is_finite() && r > 0.0),
        "density stayed positive"
    );
}

/// the analytic square-duct series solving mu lap(u) = -gpd on
/// [-a, a]^2, 40x40 odd modes (wavenumbers mm*pi/(2a), prefactor
/// 16*gpd/(mu*pi^2); the first draft of this reference had both wrong).
fn duct_u(y: f64, z: f64, gpd: f64, mu: f64) -> f64 {
    let a = 0.5f64;
    let mut s = 0.0f64;
    for m in 0..40 {
        for n in 0..40 {
            let (mm, nn) = (2 * m + 1, 2 * n + 1);
            let sign = (-1.0f64).powi(m + n);
            let am = (mm as f64) * std::f64::consts::PI / (2.0 * a);
            let bn = (nn as f64) * std::f64::consts::PI / (2.0 * a);
            s += sign / ((mm as f64) * (nn as f64) * (am * am + bn * bn))
                * (am * y).cos()
                * (bn * z).cos();
        }
    }
    s * 16.0 * gpd / (std::f64::consts::PI.powi(2) * mu)
}

#[test]
fn square_duct_profile_is_held_steady() {
    // the duct: no-slip on all four lateral faces, driven by the body
    // force matching the analytic series. initialized with the steady
    // profile and marched briefly, the profile must hold.
    let b = Bounds3d {
        xmin: 0.0,
        xmax: 0.5,
        ymin: -0.5,
        ymax: 0.5,
        zmin: -0.5,
        zmax: 0.5,
    };
    let g = grid3d(12, 16, 16, &b).unwrap();
    let mu = 0.05;
    // u_max = 0.295 * gpd * a^2 / mu: pick gpd for u_max ~ 0.05.
    let gpd = 0.0339;
    let n = g.nx * g.ny * g.nz;
    let mut u = vec![0.0f64; n];
    let mut u_max = 0.0f64;
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let c = (k * g.ny + j) * g.nx + i;
                u[c] = duct_u(g.centers_y[c], g.centers_z[c], gpd, mu);
                u_max = u_max.max(u[c]);
            }
        }
    }
    let rho = vec![1.0; n];
    let (v, w) = (vec![0.0; n], vec![0.0; n]);
    let et: Vec<f64> = u
        .iter()
        .map(|x| 1.0 / (GAMMA - 1.0) + 0.5 * x * x)
        .collect();
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &u, &v, &w, &et).unwrap();
    let mut st = ConservedState3d { rho, mx, my, mz, e };
    let walls = Walls3d {
        south: true,
        north: true,
        bottom: true,
        top: true,
        ..Default::default()
    };
    let mut t = 0.0;
    while t < 0.02 {
        let (dt, _) = advance3d_visc_rk2(
            &mut st,
            &g,
            GAMMA,
            0.5,
            true,
            ViscParams3d { mu, pr: 0.72 },
            &walls,
        )
        .unwrap();
        for mxv in st.mx.iter_mut() {
            *mxv += dt * gpd;
        }
        t += dt;
    }
    let (u_c, _, _, _) = cons_to_prim3d(&st.rho, &st.mx, &st.my, &st.mz, &st.e).unwrap();
    let i_mid = g.nx / 2;
    let mut worst = 0.0f64;
    for k in 0..g.nz {
        for j in 0..g.ny {
            let c = (k * g.ny + j) * g.nx + i_mid;
            let ana = duct_u(g.centers_y[c], g.centers_z[c], gpd, mu);
            worst = worst.max((u_c[c] - ana).abs());
        }
    }
    eprintln!(
        "duct worst profile error: {worst:.4e} ({:.2}% of u_max)",
        100.0 * worst / u_max
    );
    assert!(
        worst / u_max < 0.05,
        "duct profile must hold, worst {worst:.2e} vs u_max {u_max:.2e}"
    );
    // the cell-sampled mean/max must match the analytic field's own
    // sampled ratio (0.5699 is the integral value; center sampling
    // shifts it, so compare like against like).
    let mut sim_mean = 0.0f64;
    let mut sim_max = 0.0f64;
    for k in 0..g.nz {
        for j in 0..g.ny {
            let c = (k * g.ny + j) * g.nx + i_mid;
            sim_mean += u_c[c];
            sim_max = sim_max.max(u_c[c]);
        }
    }
    sim_mean /= (g.ny * g.nz) as f64;
    let mut ana_mean = 0.0f64;
    let mut ana_max = 0.0f64;
    for k in 0..g.nz {
        for j in 0..g.ny {
            let c = (k * g.ny + j) * g.nx + i_mid;
            let a = duct_u(g.centers_y[c], g.centers_z[c], gpd, mu);
            ana_mean += a;
            ana_max = ana_max.max(a);
        }
    }
    ana_mean /= (g.ny * g.nz) as f64;
    let sim_ratio = sim_mean / sim_max;
    let ana_ratio = ana_mean / ana_max;
    eprintln!("mean/max ratio: sim {sim_ratio:.4} vs analytic-sampled {ana_ratio:.4}");
    assert!(
        (sim_ratio - ana_ratio).abs() < 0.02,
        "duct mean/max ratio {sim_ratio} vs analytic {ana_ratio}"
    );
}

#[test]
fn viscous_operator_is_symmetric_under_y_z_swap() {
    // u = S*y between y-walls must drain exactly like u = S*z between
    // z-walls under the axis swap; guards the z-face momentum path the
    // channel test (z-invariant flow) cannot see.
    use nufor_core::{add_viscous3d, prim_to_cons3d, ViscParams3d as V3};
    let (s, mu) = (2.0f64, 0.05f64);
    let b = Bounds3d {
        xmin: 0.0,
        xmax: 0.5,
        ymin: 0.0,
        ymax: 1.0,
        zmin: 0.0,
        zmax: 1.0,
    };
    let g = grid3d(6, 8, 8, &b).unwrap();
    let n = g.nx * g.ny * g.nz;
    let rho = vec![1.0; n];
    let (v, w) = (vec![0.0; n], vec![0.0; n]);
    let idx = |i: usize, j: usize, k: usize| (k * g.ny + j) * g.nx + i;

    let mut ua = vec![0.0f64; n];
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                ua[idx(i, j, k)] = s * g.centers_y[idx(i, j, k)];
            }
        }
    }
    let et: Vec<f64> = ua
        .iter()
        .map(|x| 1.0 / (GAMMA - 1.0) + 0.5 * x * x)
        .collect();
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &ua, &v, &w, &et).unwrap();
    let mut sta = ConservedState3d {
        rho: rho.clone(),
        mx,
        my,
        mz,
        e,
    };
    let walls_y = Walls3d {
        south: true,
        north: true,
        ..Default::default()
    };
    add_viscous3d(&mut sta, &g, GAMMA, V3 { mu, pr: 0.72 }, &walls_y, 1e-3).unwrap();

    let mut ub = vec![0.0f64; n];
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                ub[idx(i, j, k)] = s * g.centers_z[idx(i, j, k)];
            }
        }
    }
    let et: Vec<f64> = ub
        .iter()
        .map(|x| 1.0 / (GAMMA - 1.0) + 0.5 * x * x)
        .collect();
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &ub, &v, &w, &et).unwrap();
    let mut stb = ConservedState3d { rho, mx, my, mz, e };
    let walls_z = Walls3d {
        bottom: true,
        top: true,
        ..Default::default()
    };
    add_viscous3d(&mut stb, &g, GAMMA, V3 { mu, pr: 0.72 }, &walls_z, 1e-3).unwrap();

    let mut worst = 0.0f64;
    for k in 0..g.nz {
        for j in 0..g.ny {
            for i in 0..g.nx {
                let d = (sta.mx[idx(i, j, k)] - stb.mx[idx(i, k, j)]).abs();
                worst = worst.max(d);
            }
        }
    }
    assert!(
        worst == 0.0,
        "y-path and z-path drains differ by {worst:.2e}"
    );
}

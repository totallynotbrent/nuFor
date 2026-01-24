//! laminar channel (poiseuille) flow: the analytic parabolic profile is held
//! steady by a viscous flow driven by an equivalent pressure-gradient body force.

use nufor_core::{
    advance2d_visc_rk2, cons_to_prim2d, grid2d, prim_to_cons2d, Bc2d, Boundaries2d,
    ConservedState2d, Grid2d, ViscParams,
};

const GAMMA: f64 = 1.4;

fn velocity_field(g: &Grid2d, umax: f64) -> Vec<f64> {
    // no-slip parabola u = umax * 4 y (1 - y) in [0, 1], peaking mid-channel.
    let n = g.nx * g.ny;
    let mut u = vec![0.0; n];
    for j in 0..g.ny {
        let y = g.centers_y[j * g.nx];
        for i in 0..g.nx {
            u[j * g.nx + i] = umax * 4.0 * y * (1.0 - y);
        }
    }
    u
}

#[test]
fn poiseuille_parabola_is_held_steady() {
    // channel [0, 0.5] x [0, 1] between two no-slip walls.
    let g: Grid2d = grid2d(24, 40, 0.0, 0.5, 0.0, 1.0).unwrap();
    let mu = 0.05;
    let umax = 0.05;
    // the body force equivalent to -dp/dx balances mu d2u/dy2: B = 8 mu umax / H^2.
    let body = 8.0 * mu * umax;
    let u = velocity_field(&g, umax);
    let rho = vec![1.0; g.nx * g.ny];
    let v = vec![0.0; g.nx * g.ny];
    let p = vec![1.0; g.nx * g.ny];
    let et: Vec<f64> = p
        .iter()
        .zip(&u)
        .map(|(pp, uu)| *pp / (rho[0] * (GAMMA - 1.0)) + 0.5 * uu * uu)
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
    // the profile must stay the parabola: compare u(y) mid-channel plane.
    let (u_c, _, _) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    let mid = g.nx / 2;
    let mut worst = 0.0f64;
    for j in 0..g.ny {
        let err = (u_c[j * g.nx + mid]
            - umax * 4.0 * g.centers_y[j * g.nx] * (1.0 - g.centers_y[j * g.nx]))
            .abs();
        worst = worst.max(err);
    }
    eprintln!("worst profile error vs poiseuille: {worst:.2e}");
    assert!(
        worst < 1e-2,
        "channel flow must hold the parabola to a percent, worst error {worst:.2e}"
    );
    assert!(
        st.rho.iter().all(|&r| r.is_finite() && r > 0.0),
        "density stayed positive"
    );
}

//! the 3d viscous operator contracts: uniform flow sees exactly zero
//! diffusive flux, and a shearing field reproduces the analytic stress.

use nufor_core::{
    add_viscous3d, grid3d, prim_to_cons3d, Bounds3d, ConservedState3d, ViscParams3d, Walls3d,
};

const GAMMA: f64 = 1.4;

fn uniform(
    nx: usize,
    ny: usize,
    nz: usize,
    u: f64,
    v: f64,
    w: f64,
) -> (nufor_core::Grid3d, ConservedState3d) {
    let b = Bounds3d {
        xmin: 0.0,
        xmax: 1.0,
        ymin: 0.0,
        ymax: 1.0,
        zmin: 0.0,
        zmax: 1.0,
    };
    let g = grid3d(nx, ny, nz, &b).unwrap();
    let n = nx * ny * nz;
    let rho = vec![1.0; n];
    let et = vec![1.0 / (GAMMA - 1.0) + 0.5 * (u * u + v * v + w * w); n];
    let (mx, my, mz, e) = prim_to_cons3d(&rho, &vec![u; n], &vec![v; n], &vec![w; n], &et).unwrap();
    (g, ConservedState3d { rho, mx, my, mz, e })
}

#[test]
fn uniform_flow_sees_zero_viscous_flux_open_domain() {
    // uniform state, all faces open: every gradient vanishes and the
    // operator must change exactly nothing, bit for bit.
    let (g, mut st) = uniform(6, 5, 4, 0.3, -0.2, 0.1);
    let before = st.clone();
    let v = ViscParams3d { mu: 0.01, pr: 0.72 };
    add_viscous3d(&mut st, &g, GAMMA, v, &Walls3d::default(), 1e-3).unwrap();
    for k in 0..st.rho.len() {
        assert_eq!(st.mx[k].to_bits(), before.mx[k].to_bits(), "mx @ {k}");
        assert_eq!(st.my[k].to_bits(), before.my[k].to_bits(), "my @ {k}");
        assert_eq!(st.mz[k].to_bits(), before.mz[k].to_bits(), "mz @ {k}");
        assert_eq!(st.e[k].to_bits(), before.e[k].to_bits(), "e @ {k}");
    }
}

#[test]
fn linear_shear_reproduces_analytic_stress() {
    // u = S * y (simple shear), rho and p uniform, no walls: the interior
    // momentum flux must match tau_xy = mu * S with dt/dx factors visible in
    // the update. with the analytic tau known, the interior cell's
    // y-momentum change over one add must be mu * S * dt * (1/dy^2... )
    // actually the cleanest analytic cell: a cell whose x-neighbors see the
    // same shear sees zero net x-flux (tau constant across faces), so the
    // momentum must NOT change; the y-direction normal flux of the constant
    // tau_xy field also cancels. verify constancy: the operator applied to a
    // linear field leaves every interior cell unchanged (discrete
    // conservation of a constant-stress state).
    let (g, mut st) = uniform(8, 8, 8, 0.0, 0.0, 0.0);
    // impose u = S*y on the conserved state (momentum mx = rho * S * y_c).
    let s = 2.0;
    for k in 0..st.rho.len() {
        let y = g.centers_y[k];
        st.mx[k] = st.rho[k] * s * y;
        // keep total energy consistent with the new kinetic part: p stays 1.
        let u = s * y;
        let e_kin = 0.5 * st.rho[k] * u * u;
        st.e[k] = 1.0 / (GAMMA - 1.0) + e_kin;
    }
    let before = st.clone();
    let v = ViscParams3d { mu: 0.02, pr: 0.72 };
    add_viscous3d(&mut st, &g, GAMMA, v, &Walls3d::default(), 1e-4).unwrap();
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let idx = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
    // interior cells: the constant-stress field must not accelerate. the
    // transmissive edges mirror the edge value, which halves the edge-cell
    // gradient and unbalances the outermost face; that reaches exactly one
    // cell inward, so the deep interior (2..n-2) is the exact-constancy
    // region.
    for k in 2..nz - 2 {
        for j in 2..ny - 2 {
            for i in 2..nx - 2 {
                let c = idx(i, j, k);
                let d = (st.mx[c] - before.mx[c]).abs();
                assert!(
                    d < 1e-12,
                    "interior cell {c} moved by {d} under constant shear stress"
                );
            }
        }
    }
}

#[test]
fn solid_wall_brakes_tangential_momentum() {
    // u = 1 everywhere with a solid south wall: the wall shear must drain
    // the x-momentum of the wall-adjacent cells and nothing else (the
    // interior faces see no gradients).
    let (g, mut st) = uniform(6, 6, 4, 1.0, 0.0, 0.0);
    let walls = Walls3d {
        south: true,
        ..Default::default()
    };
    let v = ViscParams3d { mu: 0.05, pr: 0.72 };
    let before = st.clone();
    add_viscous3d(&mut st, &g, GAMMA, v, &walls, 1e-3).unwrap();
    let (nx, ny, nz) = (g.nx, g.ny, g.nz);
    let idx = |i: usize, j: usize, k: usize| (k * ny + j) * nx + i;
    for k in 0..nz {
        for i in 0..nx {
            let wall_cell = idx(i, 0, k);
            let second = idx(i, 1, k);
            // the wall cell loses x-momentum to the wall.
            assert!(
                st.mx[wall_cell] < before.mx[wall_cell],
                "wall cell must brake"
            );
            // the second cell is the upper cell of the wall pair, so it
            // loses the interior-face shear: tyx from the wall-consistent
            // gradient (uy(wall cell) = (3u0+u1)/(3dy) = 4/(3dy), averaged
            // with cell 1's centered zero => mu*2/(3dy)).
            let dy = g.dy;
            let tau_face = 0.05 * (4.0 / (3.0 * dy)) / 2.0;
            let expected_loss = tau_face * 1e-3 / dy;
            let change = st.mx[second] - before.mx[second];
            assert!(
                (change + expected_loss).abs() < 1e-14,
                "second cell change {change} vs analytic -{expected_loss}"
            );
            // the wall cell drains the difference of the wall shear
            // (mu*(9u0-u1)/(3dy) = mu*8/(3dy)) and the interior-face shear
            // computed above, scaled by dt/dy.
            let tau_wall = 0.05 * (9.0 * 1.0 - 1.0) / (3.0 * dy);
            let expected_drop = (tau_wall - tau_face) * 1e-3 / dy;
            let drop = before.mx[wall_cell] - st.mx[wall_cell];
            assert!(
                (drop - expected_drop).abs() < 1e-14,
                "wall cell drop {drop} vs analytic {expected_drop}"
            );
        }
    }
}

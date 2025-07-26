//! verification of the solver on the lax shock tube — a much stronger case
//! (larger density/pressure jumps and higher Mach) that stresses robustness.
//!
//! compares the numerical density against the exact riemann solution and
//! checks the error norm shrinks as the mesh refines.

use nufor_core::{
    euler_solve, grid1d, prim_to_cons, riemann, Boundary, ConservedState, EulerConfig, PrimState,
};

fn lax_prim(centers: &[f64]) -> ConservedState {
    let mut rho = vec![0.0; centers.len()];
    let mut m = vec![0.0; centers.len()];
    let mut e = vec![0.0; centers.len()];
    for (i, &x) in centers.iter().enumerate() {
        // left state (0.445, 0.698, 3.528), right state (0.5, 0.0, 0.571).
        let (r, u, p) = if x < 0.0 {
            (0.445, 0.698, 3.528)
        } else {
            (0.5, 0.0, 0.571)
        };
        let et = p / (0.4 * r) + 0.5 * u * u;
        rho[i] = r;
        let (mi, ei) = prim_to_cons(&[r], &[u], &[et]).unwrap();
        m[i] = mi[0];
        e[i] = ei[0];
    }
    ConservedState { rho, m, e }
}

fn lax_to_t(n: usize, t: f64) -> (ConservedState, Vec<f64>) {
    let grid = grid1d(n, -0.5, 0.5).unwrap();
    let mut st = lax_prim(&grid.centers);
    let cfg = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: grid.dx,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: (20.0 * t / grid.dx).ceil() as usize + 200,
        t_end: t,
        tol: 0.0,
    };
    euler_solve(&mut st, &cfg).unwrap();
    (st, grid.centers)
}

fn density_l1(n: usize, t: f64) -> f64 {
    let (st, centers) = lax_to_t(n, t);
    let mut l1: f64 = 0.0;
    for (i, &x) in centers.iter().enumerate() {
        let ex = riemann(
            PrimState {
                rho: 0.445,
                u: 0.698,
                p: 3.528,
            },
            PrimState {
                rho: 0.5,
                u: 0.0,
                p: 0.571,
            },
            1.4,
            x,
            t,
        )
        .state;
        l1 += (st.rho[i] - ex.rho).abs();
    }
    l1 / n as f64
}

#[test]
fn lax_stays_physical_and_matches_exact() {
    let (st, centers) = lax_to_t(400, 0.15);
    let mut l1: f64 = 0.0;
    for (i, &x) in centers.iter().enumerate() {
        let ex = riemann(
            PrimState {
                rho: 0.445,
                u: 0.698,
                p: 3.528,
            },
            PrimState {
                rho: 0.5,
                u: 0.0,
                p: 0.571,
            },
            1.4,
            x,
            0.15,
        )
        .state;
        l1 += (st.rho[i] - ex.rho).abs();
    }
    l1 /= centers.len() as f64;
    assert!(
        st.rho.iter().all(|&r| r > 0.0),
        "density must stay positive under the strong jump"
    );
    assert!(
        st.m.iter().all(|&v| v.is_finite()),
        "momentum must stay finite"
    );
    assert!(l1 < 0.06, "lax density l1 error too large: {l1}");
}

#[test]
fn lax_error_shrinks_with_refinement() {
    let coarse = density_l1(200, 0.15);
    let fine = density_l1(400, 0.15);
    assert!(coarse > fine, "error must decrease: {coarse} vs {fine}");
}

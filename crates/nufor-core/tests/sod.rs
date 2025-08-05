//! verification of the finite-volume solver against the exact sod solution.
//!
//! runs the numerical scheme on the sod shock tube and compares the density
//! against the exact riemann solution (Toro), reporting the l1 and l-inf error
//! norms and checking that the error shrinks as the mesh refines.

use nufor_core::{
    euler_solve, grid1d, prim_to_cons, riemann, Boundary, ConservedState, EulerConfig, PrimState,
};

// the exact sod contact pressure and velocity, from the known solution.
const P_STAR_EXACT: f64 = 0.30313;
const U_STAR_EXACT: f64 = 0.92745;

fn sod_prim(centers: &[f64]) -> ConservedState {
    let mut rho = vec![0.0; centers.len()];
    let mut m = vec![0.0; centers.len()];
    let mut e = vec![0.0; centers.len()];
    for (i, &x) in centers.iter().enumerate() {
        let (r, u, et) = if x < 0.0 {
            prims(1.0, 0.0, 1.0)
        } else {
            prims(0.125, 0.0, 0.1)
        };
        rho[i] = r;
        let (mi, ei) = prim_to_cons(&[r], &[u], &[et]).unwrap();
        m[i] = mi[0];
        e[i] = ei[0];
    }
    ConservedState { rho, m, e }
}

// computes (rho, u, et) for a primitive (rho, u, p) given gamma = 1.4.
fn prims(rho: f64, u: f64, p: f64) -> (f64, f64, f64) {
    let et = p / (0.4 * rho) + 0.5 * u * u;
    (rho, u, et)
}

fn run_sod_to_t(n: usize, t: f64) -> (ConservedState, Vec<f64>) {
    let grid = grid1d(n, -0.5, 0.5).unwrap();
    let mut st = sod_prim(&grid.centers);
    let cfg = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: grid.dx,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        // a generous step budget so every mesh actually advances to time t.
        max_steps: (20.0 * t / grid.dx).ceil() as usize + 200,
        t_end: t,
        tol: 0.0,
    };
    euler_solve(&mut st, &cfg).unwrap();
    (st, grid.centers)
}

#[test]
fn exact_riemann_recovers_the_known_sod_star_state() {
    // sample well into the star region, away from the contact and the waves.
    let sol = riemann(
        PrimState {
            rho: 1.0,
            u: 0.0,
            p: 1.0,
        },
        PrimState {
            rho: 0.125,
            u: 0.0,
            p: 0.1,
        },
        1.4,
        0.15,
        0.22,
    );
    assert!(
        (sol.p_star - P_STAR_EXACT).abs() < 1e-3,
        "p* = {}",
        sol.p_star
    );
    assert!(
        (sol.u_star - U_STAR_EXACT).abs() < 1e-3,
        "u* = {}",
        sol.u_star
    );
}

#[test]
fn numerical_sod_matches_the_exact_solution() {
    // density l1 norm of the error on a moderately fine mesh.
    let (st, centers) = run_sod_to_t(400, 0.22);
    let mut l1: f64 = 0.0;
    let mut linf: f64 = 0.0;
    assert!(
        st.rho.iter().all(|&r| r > 0.0),
        "density must stay positive"
    );
    for (i, &x) in centers.iter().enumerate() {
        let ex = riemann(
            PrimState {
                rho: 1.0,
                u: 0.0,
                p: 1.0,
            },
            PrimState {
                rho: 0.125,
                u: 0.0,
                p: 0.1,
            },
            1.4,
            x,
            0.22,
        )
        .state;
        let err = (st.rho[i] - ex.rho).abs();
        l1 += err;
        linf = linf.max(err);
    }
    l1 /= centers.len() as f64;
    assert!(l1 < 0.03, "density l1 error too large: {l1}");
    assert!(linf < 0.2, "density l-inf error too large: {linf}");
}

#[test]
fn sod_error_norm_shrinks_with_mesh_refinement() {
    let coarse = error_l1(200);
    let fine = error_l1(400);
    assert!(
        coarse > fine,
        "error must decrease with refinement: {coarse} > {fine}"
    );
}

fn error_l1(n: usize) -> f64 {
    let (st, centers) = run_sod_to_t(n, 0.22);
    let mut l1: f64 = 0.0;
    for (i, &x) in centers.iter().enumerate() {
        let ex = riemann(
            PrimState {
                rho: 1.0,
                u: 0.0,
                p: 1.0,
            },
            PrimState {
                rho: 0.125,
                u: 0.0,
                p: 0.1,
            },
            1.4,
            x,
            0.22,
        )
        .state;
        l1 += (st.rho[i] - ex.rho).abs();
    }
    l1 / n as f64
}

//! the verification ladder: constant state, uniform advection, a stationary
//! shock, and an isentropic expansion, each exercised through the solver.
//!
//! the ladder is the project's structured guarantee that the scheme is not
//! just stable but physically faithful one rung at a time.

use nufor_core::{
    eos_pressure, euler_solve, riemann, Boundary, ConservedState, EulerConfig, PrimState,
};

fn cfg(dx: f64, left: Boundary, right: Boundary, max_steps: usize, t_end: f64) -> EulerConfig {
    EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx,
        left,
        right,
        max_steps,
        t_end,
        tol: 0.0,
    }
}

// uniform moving state: every cell the same, so every face flux is the same
// and the state must not change regardless of how far it is advanced.
#[test]
fn uniform_state_advects_without_change() {
    let n = 64;
    let dx = 1.0 / n as f64;
    let (m0, e0) = nufor_core::prim_to_cons(&[1.0], &[2.0], &[5.0]).unwrap();
    let mut st = ConservedState {
        rho: vec![1.0; n],
        m: vec![m0[0]; n],
        e: vec![e0[0]; n],
    };
    let c = cfg(dx, Boundary::Transmissive, Boundary::Transmissive, 200, 1.0);
    let res = euler_solve(&mut st, &c).unwrap();
    assert!(res.steps > 0);
    let max_delta = st
        .rho
        .iter()
        .zip(&st.m)
        .zip(&st.e)
        .map(|((&r, &m), &e)| {
            (r - 1.0)
                .abs()
                .max((m - m0[0]).abs())
                .max((e - e0[0]).abs())
        })
        .fold(0.0_f64, f64::max);
    assert!(
        max_delta < 1e-9,
        "uniform advection must not distort the state, max delta {max_delta}"
    );
}

// a stationary normal shock built from the rankine-hugoniot relations: the two
// plateaus either side must persist through the run.
#[test]
fn stationary_shock_keeps_its_plateaus() {
    // upstream rho=1, p=1, mach 2.0 shock at rest in this frame.
    let gamma: f64 = 1.4;
    let m1: f64 = 2.0;
    let (r1, p1): (f64, f64) = (1.0, 1.0);
    let a1 = (gamma * p1 / r1).sqrt();
    let u1 = m1 * a1;
    // post-shock state from the normal-shock relations.
    let r4 = r1 * ((gamma + 1.0) * m1 * m1) / ((gamma - 1.0) * m1 * m1 + 2.0);
    let p4 = p1 * (2.0 * gamma * m1 * m1 - (gamma - 1.0)) / (gamma + 1.0);
    let u4 = u1 * r1 / r4;
    let n = 200;
    let dx = 1.0 / n as f64;
    let mut rho = Vec::with_capacity(n);
    let mut m = Vec::with_capacity(n);
    let mut e = Vec::with_capacity(n);
    for i in 0..n {
        let (r, u, p) = if i < n / 2 {
            (r1, u1, p1)
        } else {
            (r4, u4, p4)
        };
        let et = p / ((gamma - 1.0) * r) + 0.5 * u * u;
        rho.push(r);
        let (mi, ei) = nufor_core::prim_to_cons(&[r], &[u], &[et]).unwrap();
        m.push(mi[0]);
        e.push(ei[0]);
    }
    let mut st = ConservedState { rho, m, e };
    let c = cfg(dx, Boundary::Transmissive, Boundary::Transmissive, 400, 1.0);
    euler_solve(&mut st, &c).unwrap();
    assert!(st.rho.iter().all(|&r| r > 0.0));
    // plateaus far from the shock survive.
    let left_p = (0..8).map(|i| st.rho[i]).sum::<f64>() / 8.0;
    let right_p = (n - 8..n).map(|i| st.rho[i]).sum::<f64>() / 8.0;
    assert!(
        (left_p - r1).abs() / r1 < 0.02,
        "left plateau drifted: {left_p}"
    );
    assert!(
        (right_p - r4).abs() / r4 < 0.02,
        "right plateau drifted: {right_p}"
    );
    // the jump is still present.
    let max_rho = st.rho.iter().cloned().fold(0.0_f64, f64::max);
    assert!(max_rho - st.rho.iter().cloned().fold(f64::INFINITY, f64::min) > r4 * 0.5);
}

// an isentropic expansion (the left rarefaction of the sod case): the quantity
// p * rho^-gamma must be constant through the fan.
#[test]
fn isentropic_expansion_keeps_a_constant_p_rho_gamma() {
    let gamma = 1.4;
    let n = 400;
    let dx = 1.0 / n as f64;
    let mut st = ConservedState {
        rho: vec![0.0; n],
        m: vec![0.0; n],
        e: vec![0.0; n],
    };
    for i in 0..n {
        let (r, u, p) = if (i as f64 + 0.5) * dx < 0.5 {
            (1.0, 0.0, 1.0)
        } else {
            (0.125, 0.0, 0.1)
        };
        st.rho[i] = r;
        let et = p / ((gamma - 1.0) * r);
        let (mi, ei) = nufor_core::prim_to_cons(&[r], &[u], &[et]).unwrap();
        st.m[i] = mi[0];
        st.e[i] = ei[0];
    }
    let c = cfg(
        dx,
        Boundary::Transmissive,
        Boundary::Transmissive,
        2000,
        0.2,
    );
    euler_solve(&mut st, &c).unwrap();
    // find the interior cells that sit in the rarefaction (between the leading
    // edge and the contact, well away from the boundaries) and check the invariant.
    let u = |i: usize| st.m[i] / st.rho[i];
    let et = |i: usize| st.e[i] / st.rho[i];
    let p = |i: usize| eos_pressure(gamma, &[st.rho[i]], &[et(i)], &[u(i)]).unwrap()[0];
    let mut vals = Vec::new();
    for i in 20..n / 2 - 5 {
        let x = (i as f64 + 0.5) * dx;
        // only cells left of the contact (x < ~0.5+0.927*... ) and inside the fan.
        if x < 0.5 + 0.5 * 0.2 && x > 0.5 - 0.5 * 2.0 * 0.2 {
            vals.push(p(i) * st.rho[i].powf(-gamma));
        }
    }
    assert!(!vals.is_empty());
    let mean = vals.iter().sum::<f64>() / vals.len() as f64;
    // the sod left state p*rho^-g = 1.0; the fan must stay within a few percent.
    assert!(
        (mean - 1.0).abs() < 0.05,
        "isentropic invariant drifted: {mean}"
    );
    let spread = vals
        .iter()
        .map(|&v| (v - mean).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        spread < 0.03,
        "isentropic invariant not constant: max spread {spread}"
    );
}

// the constant-state rung: a reference against the exact flat solution.
#[test]
fn constant_state_matches_exact_flat_solution() {
    let sol = riemann(
        PrimState {
            rho: 1.0,
            u: 0.0,
            p: 1.0,
        },
        PrimState {
            rho: 1.0,
            u: 0.0,
            p: 1.0,
        },
        1.4,
        0.3,
        0.5,
    );
    assert!((sol.state.rho - 1.0).abs() < 1e-12);
    assert!(sol.state.u.abs() < 1e-12);
    assert!((sol.state.p - 1.0).abs() < 1e-12);
}

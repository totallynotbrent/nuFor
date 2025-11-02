//! diagnostics: physical-validity scanning and solver termination reasons.

use nufor_core::{
    advance, check_physical, euler_solve, prim_to_cons, Boundary, ConservedState, EulerConfig,
    TerminationReason,
};

fn sod(n: usize) -> ConservedState {
    let mut s = ConservedState {
        rho: vec![0.0; n],
        m: vec![0.0; n],
        e: vec![0.0; n],
    };
    for i in 0..n {
        let cx = (i as f64 + 0.5) / n as f64;
        let (r, u, p) = if cx < 0.5 {
            (1.0, 0.0, 1.0)
        } else {
            (0.125, 0.0, 0.1)
        };
        s.rho[i] = r;
        let et = p / (0.4 * r) + 0.5 * u * u;
        let (m, e) = prim_to_cons(&[r], &[u], &[et]).unwrap();
        s.m[i] = m[0];
        s.e[i] = e[0];
    }
    s
}

#[test]
fn scan_flags_non_finite_density() {
    let mut s = sod(16);
    s.rho[7] = f64::NAN;
    let chk = check_physical(&s, 1.4);
    assert!(!chk.ok);
    assert_eq!(chk.bad_cell, Some(7));
}

#[test]
fn scan_flags_negative_density() {
    let mut s = sod(16);
    s.rho[3] = -0.5;
    let chk = check_physical(&s, 1.4);
    assert!(!chk.ok);
    assert_eq!(chk.bad_cell, Some(3));
}

#[test]
fn scan_flags_negative_pressure() {
    // a cell moving so fast that its kinetic energy exceeds the total energy.
    let mut s = sod(16);
    let (r, u, et) = (1.0, 10.0, 45.0); // 0.5 u^2 = 50 > 45 -> p < 0
    let (m, e) = prim_to_cons(&[r], &[u], &[et]).unwrap();
    s.m[9] = m[0];
    s.e[9] = e[0];
    let chk = check_physical(&s, 1.4);
    assert!(!chk.ok);
    assert_eq!(chk.bad_cell, Some(9));
    assert!(chk.min_p < 0.0);
}

#[test]
fn scan_passes_a_healthy_sod() {
    let s = sod(64);
    let chk = check_physical(&s, 1.4);
    assert!(chk.ok);
    assert_eq!(chk.bad_cell, None);
    assert!((chk.min_rho - 0.125).abs() < 1e-12);
    assert!(chk.min_p > 0.0);
}

#[test]
fn reason_is_converged_when_tolerance_is_met() {
    let mut s = sod(50);
    let cfg = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: 1.0 / 50.0,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: 10_000,
        t_end: f64::INFINITY,
        tol: 1.0, // the first-step residual is well below this
    };
    let r = euler_solve(&mut s, &cfg).unwrap();
    assert_eq!(r.reason, TerminationReason::Converged);
    assert!(r.converged);
    assert!(r.steps >= 1);
}

#[test]
fn reason_is_time_end_when_t_is_short() {
    let mut s = sod(50);
    let cfg = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: 1.0 / 50.0,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: 10_000,
        t_end: 1e-6,
        tol: 0.0,
    };
    let r = euler_solve(&mut s, &cfg).unwrap();
    assert_eq!(r.reason, TerminationReason::TimeEnd);
    assert!(r.time >= 1e-6);
}

#[test]
fn reason_is_max_steps_when_budget_runs_out() {
    let mut s = sod(50);
    let cfg = EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: 1.0 / 50.0,
        left: Boundary::Transmissive,
        right: Boundary::Transmissive,
        max_steps: 1,
        t_end: f64::INFINITY,
        tol: 0.0,
    };
    let r = euler_solve(&mut s, &cfg).unwrap();
    assert_eq!(r.reason, TerminationReason::MaxSteps);
    assert_eq!(r.steps, 1);
}

#[test]
fn scan_catches_a_run_that_goes_non_physical() {
    // a single-cell velocity spike under open boundaries; the escaping stream
    // drives that cell toward vacuum and the density turns negative.
    let n = 50;
    let mut s = sod(n);
    let mid = n / 2;
    let (r, u, p) = (1.0, 8.0, 0.02);
    let et = p / (0.4 * r) + 0.5 * u * u;
    let (m, e) = prim_to_cons(&[r], &[u], &[et]).unwrap();
    s.m[mid] = m[0];
    s.e[mid] = e[0];
    let mut blew = false;
    for _ in 0..4000 {
        if advance(
            &mut s,
            1.4,
            0.5,
            1.0 / n as f64,
            Boundary::Transmissive,
            Boundary::Transmissive,
        )
        .is_err()
        {
            blew = true;
            break;
        }
        let chk = check_physical(&s, 1.4);
        if !chk.ok {
            blew = true;
            assert!(chk.bad_cell.is_some());
            break;
        }
    }
    assert!(blew, "the spike should eventually go non-physical");
}

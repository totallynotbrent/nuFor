//! integration tests for the 1d euler finite-volume time-marching solver.
//!
//! the solver reuses the fortran kernels for cfl, conversion, and hll flux and
//! adds the conservative update, ghost-cell boundaries, residual, and logging in rust.

use nufor_core::{
    advance, cfl_dt, euler_solve, prim_to_cons, Boundary, ConservedState, Error, EulerConfig,
};

// sod shock-tube primitives: left (1,0,2.5), right (0.125,0,2.0) at gamma=1.4.
fn sod_state(left: usize, right: usize) -> ConservedState {
    let mut rho = vec![0.0; left + right];
    let mut m = vec![0.0; left + right];
    let mut e = vec![0.0; left + right];
    for i in 0..left {
        let (mi, ei) = prim_to_cons(&[1.0], &[0.0], &[2.5]).unwrap();
        m[i] = mi[0];
        e[i] = ei[0];
        rho[i] = 1.0;
    }
    for i in left..left + right {
        let (mi, ei) = prim_to_cons(&[0.125], &[0.0], &[2.0]).unwrap();
        m[i] = mi[0];
        e[i] = ei[0];
        rho[i] = 0.125;
    }
    ConservedState { rho, m, e }
}

fn total_mass(state: &ConservedState, dx: f64) -> f64 {
    state.rho.iter().map(|&r| r * dx).sum()
}

fn total_energy(state: &ConservedState, dx: f64) -> f64 {
    state.e.iter().map(|&v| v * dx).sum()
}

// reusable sod config; callers override the boundary and stopping rules.
fn sod_cfg(left: Boundary, right: Boundary, max_steps: usize, t_end: f64) -> EulerConfig {
    EulerConfig {
        gamma: 1.4,
        cfl: 0.5,
        dx: 1.0,
        left,
        right,
        max_steps,
        t_end,
        tol: 1e-6,
    }
}

// gamma=2, uniform (rho=1,u=0,et=8): the hll face flux is constant so no change.
fn uniform_state() -> ConservedState {
    let (m, e) = prim_to_cons(&[1.0], &[0.0], &[8.0]).unwrap();
    ConservedState {
        rho: vec![1.0; 8],
        m: vec![m[0]; 8],
        e: vec![e[0]; 8],
    }
}

#[test]
fn uniform_state_stays_at_rest() {
    let mut st = uniform_state();
    let (_, resid) = advance(
        &mut st,
        2.0,
        0.5,
        1.0,
        Boundary::Transmissive,
        Boundary::Transmissive,
    )
    .unwrap();
    assert!(
        resid < 1e-12,
        "uniform state must have ~zero residual, got {resid}"
    );
    for i in 0..st.rho.len() {
        assert!((st.rho[i] - 1.0).abs() < 1e-12);
        assert!(st.m[i].abs() < 1e-12);
    }
}

#[test]
fn reflective_walls_conserve_mass_and_energy() {
    let mut st = sod_state(100, 100);
    let dx = 1.0;
    let m0 = total_mass(&st, dx);
    let e0 = total_energy(&st, dx);
    let cfg = sod_cfg(Boundary::Reflective, Boundary::Reflective, 400, 10.0);
    let res = euler_solve(&mut st, &cfg).unwrap();
    assert!(res.steps > 0);
    // a closed tube conserves mass and total energy to rounding.
    let dm = (total_mass(&st, dx) - m0).abs() / m0;
    let de = (total_energy(&st, dx) - e0).abs() / e0;
    assert!(dm < 1e-10, "mass not conserved: rel {dm}");
    assert!(de < 1e-10, "energy not conserved: rel {de}");
    // the flow stays physical.
    assert!(st.rho.iter().all(|&r| r > 0.0));
}

#[test]
fn sod_evolves_physically_on_transmissive_ends() {
    let mut st = sod_state(100, 100);
    let cfg = sod_cfg(Boundary::Transmissive, Boundary::Transmissive, 400, 5.0);
    let res = euler_solve(&mut st, &cfg).unwrap();
    assert!(res.steps > 0);
    assert!(res.time > 0.0);
    assert!(
        st.rho.iter().all(|&r| r > 0.0),
        "density must stay positive"
    );
    // expansion should have carried density low but finite, and a shock region should appear.
    let min = st.rho.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_m2 =
        st.m.iter()
            .map(|m| (m / 1.0).powi(2))
            .fold(0.0_f64, f64::max);
    assert!(
        min > 0.0 && min < 1.0,
        "sod expansion must lower density, got min {min}"
    );
    assert!(max_m2.is_finite());
}

#[test]
fn first_step_uses_the_cfl_time_step() {
    let mut st = sod_state(50, 50);
    // cfl on the initial state is what advance should use for its first step.
    let step = cfl_dt(1.4, 0.5, 1.0, &st.rho, &st.m, &st.e).unwrap();
    let (dt, _) = advance(
        &mut st,
        1.4,
        0.5,
        1.0,
        Boundary::Transmissive,
        Boundary::Transmissive,
    )
    .unwrap();
    assert!((dt - step.dt).abs() < 1e-12);
    assert!(dt > 0.0);
}

#[test]
fn single_cell_is_rejected() {
    let mut st = ConservedState {
        rho: vec![1.0; 1],
        m: vec![0.0; 1],
        e: vec![1.0; 1],
    };
    assert!(matches!(
        advance(
            &mut st,
            1.4,
            0.5,
            1.0,
            Boundary::Transmissive,
            Boundary::Transmissive
        ),
        Err(Error::InvalidArgs)
    ));
}

#[test]
fn log_time_advances_by_dt() {
    let mut st = sod_state(50, 50);
    let cfg = sod_cfg(Boundary::Reflective, Boundary::Reflective, 100, 10.0);
    let res = euler_solve(&mut st, &cfg).unwrap();
    let mut t = 0.0;
    for row in &res.log {
        assert!((row.time - t - row.dt).abs() < 1e-12);
        t = row.time;
    }
    assert!((t - res.time).abs() < 1e-12);
}

//! park two-temperature rung gates.
//!
//! 1. landau-teller box, frozen translational temperature: the relaxed
//!    e_v must track the analytic curve e_eq(T)[1-exp(-t/tau)] to
//!    tight tolerance.
//! 2. coupled box (translational temperature drifting as the
//!    vibrational pool fills): total energy per cell must be conserved
//!    exactly.
//! 3. t_v inversion: t_vib_of_e_v(e_vib_of_t_v(t)) == t roundtrip.
use nufor_core::grid2d;
use nufor_core::park2t::{e_vib_of_t_v, park2t_substep, t_vib_of_e_v, tau_vt_n2, Park2tState};
use nufor_core::ConservedState2d;

fn box_state(g: &nufor_core::Grid2d, rho: f64, t_tr: f64, u: f64) -> ConservedState2d {
    let n = g.nx * g.ny;
    let ein = t_tr / (1.4 - 1.0) * 8.314462 / 28.016e-3;
    let et = ein + 0.5 * u * u;
    ConservedState2d {
        rho: vec![rho; n],
        mx: vec![rho * u; n],
        my: vec![0.0; n],
        e: vec![rho * et; n],
    }
}

#[test]
fn park2t_tv_inversion_roundtrip() {
    for tv in [300.0, 1000.0, 3393.0, 8000.0, 30000.0] {
        let ev = e_vib_of_t_v(tv);
        let back = t_vib_of_e_v(ev);
        let rel = (back - tv).abs() / tv;
        assert!(rel < 1e-9, "tv {tv} roundtrip rel {rel:.2e}");
    }
    // zero energy reads as the cold floor
    assert_eq!(t_vib_of_e_v(0.0), 1.0);
}

#[test]
fn park2t_landau_teller_frozen_t() {
    // frozen-t box: n2 at rho = 1e-3, T = 6000 K, e_v(0) = the 300 K
    // value. the analytic relaxed curve with tau from the millikan-white
    // + cutoff model: e(t) = e_eq + (e0 - e_eq) exp(-t/tau).
    let g = grid2d(8, 8, 0.0, 1e-6, 0.0, 1e-6).unwrap();
    let n = g.nx * g.ny;
    let solid = vec![false; n];
    let rho = 1.0e-3;
    let t_tr = 6000.0;
    let e0 = e_vib_of_t_v(300.0);
    let mut park = Park2tState::from_uniform(&g, e0);
    let p = (1.4 - 1.0) * rho * (t_tr / 0.4 * 8.314462 / 28.016e-3);
    let n_den = rho / 28.016e-3 * 6.022e23;
    let tau = tau_vt_n2(t_tr, p, n_den);

    // march in 200 substeps to t = 3 tau, holding T frozen by putting
    // the draw back (the frozen-t experiment): run the real substep and
    // compare against the ANALYTIC curve WITH the draw: the frozen-t
    // gate holds only if the draw is negligible; at these conditions
    // e_eq(6000) vs the translational pool: check the draw fraction.
    let e_eq = e_vib_of_t_v(t_tr);
    // the thermostat refunds the draw after each substep, so the
    // analytic constant-T curve holds as long as the per-substep draw
    // is small relative to the pool.
    let ein_tr = t_tr / 0.4 * 8.314462 / 28.016e-3;
    let dt = tau * 3.0 / 600.0;
    let per_step_draw = (e_eq - e0) * (1.0 - (-dt / tau).exp()) / ein_tr;
    assert!(
        per_step_draw < 5e-3,
        "frozen-t precondition: per-substep draw {per_step_draw:.3e}"
    );

    let ein0_tr = t_tr / 0.4 * 8.314462 / 28.016e-3;
    for _ in 0..600 {
        // frozen-t protocol: reset the translational pool to exactly
        // t_tr each substep (a thermostat); the analytic curve assumes
        // constant T and constant p in tau. the conserved total picks
        // up the injected energy by design of the protocol.
        let mut st = box_state(&g, rho, t_tr, 0.0);
        for k in 0..n {
            st.e[k] = rho * (ein0_tr + park.e_v[k]);
        }
        park2t_substep(&mut st, &mut park, &g, 1.4, dt, &solid, false).unwrap();
    }
    let t_total = tau * 3.0;
    let e_analytic = e_eq + (e0 - e_eq) * (-t_total / tau).exp();
    for k in 0..n {
        let rel = (park.e_v[k] - e_analytic).abs() / e_analytic.abs().max(1e-30);
        assert!(rel < 5e-3, "cell {k}: e_v rel drift {rel:.3e}");
    }
}

#[test]
fn park2t_coupled_box_energy_conserved() {
    // no thermostat: the vibrational fill draws from translation. total
    // energy per cell (e_trans + e_vib) must be invariant to roundoff.
    let g = grid2d(8, 8, 0.0, 1e-6, 0.0, 1e-6).unwrap();
    let n = g.nx * g.ny;
    let solid = vec![false; n];
    let rho = 1.0e-3;
    let mut st = box_state(&g, rho, 6000.0, 0.0);
    let e0 = e_vib_of_t_v(300.0);
    let mut park = Park2tState::from_uniform(&g, e0);

    let ein0 = st.e[0] / rho;
    let p = (1.4 - 1.0) * rho * (t_of_ein(ein0));
    let n_den = rho / 28.016e-3 * 6.022e23;
    let tau = tau_vt_n2(6000.0, p, n_den);
    let dt = tau * 0.5;
    for _ in 0..200 {
        park2t_substep(&mut st, &mut park, &g, 1.4, dt, &solid, false).unwrap();
    }
    // the coupled equilibrium: e_v = e_vib(T), T = (ein0 - e_v)(g-1)/R.
    let r_n2 = 8.314462 / 28.016e-3;
    let (mut lo, mut hi) = (0.0f64, ein0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let t_eq = (ein0 - mid) * 0.4 / r_n2;
        let x = 3393.0 / t_eq.max(1.0);
        let ev = if x < 40.0 {
            r_n2 * 3393.0 / ((x.exp() - 1.0).max(1e-12))
        } else {
            0.0
        };
        if ev > mid {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let e_star = 0.5 * (lo + hi);
    for k in 0..n {
        let u = st.mx[k] / st.rho[k];
        let ein = st.e[k] / st.rho[k] - 0.5 * u * u;
        // the conserved total: ein (translational + vibrational inside
        // state.e) plus e_v ONCE — state.e already carries the full pool,
        // so the invariant is ein itself.
        let rel_e = (ein - ein0).abs() / ein0;
        assert!(rel_e < 1e-12, "cell {k}: conserved ein drift {rel_e:.3e}");
        let rel_v = (park.e_v[k] - e_star).abs() / e_star;
        assert!(
            rel_v < 1e-3,
            "cell {k}: e_v vs coupled fixed point {rel_v:.3e}"
        );
    }
}

fn t_of_ein(ein: f64) -> f64 {
    ein * 0.4 / (8.314462 / 28.016e-3)
}

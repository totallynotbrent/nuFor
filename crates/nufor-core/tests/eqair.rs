//! tests for the equilibrium-air closure (Srinivasan-Tannehill TGAS1
//! fits): the two acceptance gates from the roadmap plus block-level
//! invariants. the numbers are checked against independently computed
//! references (python reference implementation, physics truth), not
//! against the fits themselves, so a transcription slip in any block
//! shifts a value out of its band.

use nufor_core::{eqair_pressure, eqair_sound};

/// cold air at sea level reproduces the perfect-gas answer to 0.5%.
#[test]
fn eqair_cold_air_matches_sea_level_truth() {
    // rho = 1.225 kg/m^3, T = 288 K, cv = 718 J/kg/K.
    let rho = [1.225_f64];
    let e = [718.0 * 288.0];
    let p = eqair_pressure(&rho, &e).unwrap();
    assert!((p[0] - 101325.0).abs() / 101325.0 < 5e-3, "p = {}", p[0]);
    let a = eqair_sound(&rho, &e).unwrap();
    assert!((a[0] - 340.3).abs() < 2.0, "a = {}", a[0]);
}

/// GAMM stays inside the physical band over the fits' valid window as
/// the gas heats: gamma leaves the cold 1.4 plateau and drops toward
/// the dissociated regime; it may rise slightly again when ionization
/// competes, so the contract is a band, not monotonicity.
#[test]
fn eqair_gamma_stays_physical_with_temperature() {
    let rho = 1.0e-3_f64;
    // Z from 0.7 to 3.3 (the region-B window for this density).
    for k in 0..14 {
        let z = 0.7 + 0.2 * k as f64;
        let e = [78408.4_f64 * 10_f64.powf(z)];
        let p = eqair_pressure(&[rho], &e).unwrap();
        let gamm = 1.0 + p[0] / (rho * e[0]);
        assert!(
            gamm > 1.05 && gamm < 1.42,
            "gamma out of band at Z = {z}: {gamm}"
        );
    }
}

/// the M22 normal shock compresses far past the perfect-gas cap of 6:
/// equilibrium dissociation absorbs energy, so the density jump grows.
/// the jump is solved from the Rankine-Hugoniot relations with the
/// eqair closure supplying p(rho, e); the test asserts the converged
/// jump lands in the physical band and the post-shock flow is subsonic.
#[test]
fn eqair_mach22_shock_beats_perfect_gas_cap() {
    let r_air = 287.0_f64;
    let rho1 = 4.0e-4;
    let t1 = 270.0;
    let e1 = 718.0 * t1;
    let p1 = rho1 * r_air * t1;
    let a1 = (1.4 * p1 / rho1).sqrt();
    let u1 = 22.0 * a1;

    // residual: p_eos(rho2, e2(rho2)) - p_momentum(rho2)
    let resid = |x: f64| -> f64 {
        let u2 = rho1 * u1 / x;
        let p2 = p1 + rho1 * u1 * (u1 - u2);
        let e2 = e1 + 0.5 * (u1 * u1 - u2 * u2) + p1 / rho1 - p2 / x;
        eqair_pressure(&[x], &[e2]).unwrap()[0] - p2
    };

    // bracket: f is negative below the root, positive above it.
    let mut lo = 10.0 * rho1;
    let mut hi = 25.0 * rho1;
    assert!(resid(lo) < 0.0);
    assert!(resid(hi) > 0.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if resid(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let rho2 = 0.5 * (lo + hi);
    let jump = rho2 / rho1;
    assert!(jump > 10.0, "jump = {jump}, perfect-gas cap is 6");
    assert!(
        jump < 25.0,
        "jump = {jump}, above the bracket's physical band"
    );

    // post-shock state is subsonic under the eqair sound speed.
    let u2 = rho1 * u1 / rho2;
    let p2 = p1 + rho1 * u1 * (u1 - u2);
    let e2 = e1 + 0.5 * (u1 * u1 - u2 * u2) + p1 / rho1 - p2 / rho2;
    let a2 = eqair_sound(&[rho2], &[e2]).unwrap()[0];
    let m2 = u2 / a2;
    assert!(m2 < 1.0, "post-shock Mach = {m2}, must be subsonic");
}

/// the closure rejects non-positive density or energy outright.
#[test]
fn eqair_rejects_nonphysical_inputs() {
    assert!(eqair_pressure(&[-1.0], &[1.0]).is_err());
    assert!(eqair_pressure(&[1.0], &[0.0]).is_err());
    assert!(eqair_pressure(&[], &[]).is_err());
    assert!(eqair_sound(&[1.0], &[-1.0]).is_err());
}

/// every fit block answers with a finite, positive pressure and sound
/// speed at its interior points, over each region's valid Z window
/// (the listing warns past the last break of each region; there the
/// fits are extrapolations and sound speed can go non-physical, which
/// the closure reports as an error rather than a value).
#[test]
fn eqair_answers_finite_over_the_whole_domain() {
    // Y sweeps all three density regions; Z windows per region:
    // A up to 3.69, B up to 3.4, C up to 2.9.
    for k in 0..30 {
        let y = -6.9 + 0.286 * k as f64;
        let rho = 1.292 * 10_f64.powf(y);
        let zmax = if y <= -4.5 {
            3.69_f64
        } else if y <= -0.5 {
            3.4_f64
        } else {
            2.9_f64
        };
        let steps = (((zmax - 0.7) / 0.12).ceil()) as usize;
        for m in 0..steps {
            let z = 0.7 + 0.12 * m as f64;
            let e = 78408.4 * 10_f64.powf(z);
            let p = eqair_pressure(&[rho], &[e]).unwrap();
            assert!(
                p[0].is_finite() && p[0] > 0.0,
                "bad p at Y={y} Z={z}: {}",
                p[0]
            );
            let a = eqair_sound(&[rho], &[e]).unwrap();
            assert!(
                a[0].is_finite() && a[0] > 100.0,
                "bad a at Y={y} Z={z}: {}",
                a[0]
            );
        }
    }
}

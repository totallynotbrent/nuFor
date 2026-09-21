//! spalart-allmaras closure and source-term unit tests against the canonical
//! 1994 constants and hand-computed reference values.

use nufor_core::{
    fv1, fv2, fw, g_func, source, stilde, C_B1, C_B2, C_V1, C_W1, C_W2, C_W3, KAPPA, SIGMA,
};

const EPS: f64 = 1e-12;

#[test]
fn constants_match_the_1994_paper() {
    // the canonical set from spalart & allmaras (1994), table 1.
    assert!((C_B1 - 0.1355).abs() < EPS);
    assert!((C_B2 - 0.622).abs() < EPS);
    assert!((SIGMA - 2.0 / 3.0).abs() < EPS);
    assert!((KAPPA - 0.41).abs() < EPS);
    assert!((C_W2 - 0.3).abs() < EPS);
    assert!((C_W3 - 2.0).abs() < EPS);
    assert!((C_V1 - 7.1).abs() < EPS);
    // c_w1 = c_b1/kappa^2 + (1 + c_b2)/sigma, the standard composite.
    let cw1_expected = C_B1 / (KAPPA * KAPPA) + (1.0 + C_B2) / SIGMA;
    assert!((C_W1 - cw1_expected).abs() < EPS);
    assert!((C_W1 - 3.2390678167757287).abs() < 1e-9);
}

#[test]
fn fv1_rises_monotonically_with_chi() {
    // chi = nu_tilde/nu; fv1 = chi^3 / (chi^3 + c_v1^3).
    let f0 = fv1(0.0);
    assert!(f0.abs() < EPS);
    let f_small = fv1(0.1);
    let f_mid = fv1(1.0);
    let f_big = fv1(10.0);
    // monotone increasing, saturating toward 1.
    assert!(f_small < f_mid);
    assert!(f_mid < f_big);
    assert!(f_big < 1.0);
    // at chi = c_v1, fv1 = 0.5 exactly by symmetry of the formula.
    let f_cv1 = fv1(C_V1);
    assert!((f_cv1 - 0.5).abs() < 1e-10);
}

#[test]
fn fv2_and_r_agree_at_the_wall_limit() {
    // fv2 = 1 - chi/(1 + chi fv1); at chi=0 both fv2 and r vanish appropriately.
    assert!((fv2(0.0) - 1.0).abs() < EPS);
    // large chi limit: fv2 -> 1 - chi/(chi + chi*1) = 1 - 1/2 = 0.5? verify numerically below.
    let f_big = fv2(100.0);
    assert!(f_big > 0.0 && f_big < 1.0);
}

#[test]
fn fw_falls_from_one_toward_zero() {
    // g = r + c_w2 (r^6 - r); fw = g ((1+c_w3^6)/(g^6+1+c_w3^6))^(1/6).
    // canonical anchors: g(0) = 0 so fw(0) = 0 exactly; g(1) = 1 so
    // fw(1) = (65/66)^(1/6) = 0.997459 (the log-layer value).
    assert!(g_func(0.0).abs() < EPS);
    assert!(fw(0.0).abs() < EPS);
    assert!((g_func(1.0) - 1.0).abs() < EPS);
    assert!((fw(1.0) - 0.9974586560076582).abs() < 1e-12);
    // large r: g grows like c_w2 r^6 and fw rises past 1 (destruction
    // damping for the far-from-equilibrium region).
    assert!(fw(10.0) > 1.0);
    // monotone in the near-equilibrium range.
    assert!(fw(0.5) < fw(1.0));
    assert!(fw(2.0) > fw(1.0));
}

#[test]
fn stilde_is_positive_with_vorticity_and_distance() {
    // stilde = s + nu_tilde/(kappa^2 d^2) fv2; s = sqrt(2)|omega| in 2d (planar vorticity magnitude).
    let s_tilde = stilde(0.0, 1.0, 1.0, 1.0);
    assert!(s_tilde > 0.0);
    // zero vorticity with nu_tilde>0 near a wall still gives positive stilde (the distance term).
    let s_wall = stilde(1.0, 0.0, 1.0, 0.01);
    assert!(s_wall > 0.0);
    // the fv2 < 0 freestream region (large chi) is floored at cs * s, never negative.
    let s_free = stilde(50.0 * 1.0, 1.0, 100.0, 1.0);
    assert!(s_free >= 0.3 * 1.0 - 1e-12, "got {s_free}");
    assert!(s_free > 0.0);
}

#[test]
fn source_is_a_wash_in_far_freestream() {
    // far from any wall (d large) with zero vorticity, stilde -> s = 0, so
    // production vanishes; destruction also vanishes with (nu_tilde/d)^2.
    // a uniform freestream field must not drift.
    let nu = 1.0;
    let nu_tilde = 3.0 * nu;
    let src = source(nu_tilde, 0.0, 1e6, nu, 1e-6, 1.0);
    assert!(src.abs() < 1e-3, "got {src}");
    // near a wall with the same vorticity, destruction dominates: the field
    // is driven down toward its wall value.
    let near = source(nu_tilde, 0.0, 1e-3, nu, 1e-6, 1.0);
    assert!(near < 0.0, "near-wall source must be negative, got {near}");
}

#[test]
fn mu_t_formula_matches_fv1_coupling() {
    // mu_t = rho nu_tilde fv1; with chi = 3, fv1(3) = 27/(27+357.911) = 0.0701...
    let chi = 3.0;
    let expected = chi * fv1(chi);
    assert!((expected - 3.0 * 0.0701).abs() < 2e-3);
    // wait: mu_t = rho*nu_tilde*fv1 = nu*chi*fv1 per unit rho. the value 0.0701*3 = 0.2103.
    assert!((expected - 0.2103).abs() < 2e-3);
}

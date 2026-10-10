//! the sphere-cone geometry contract: the reference capsule meridian must
//! reproduce the trade-study numbers exactly, and the signed-distance
//! body must be consistent with them (negative inside, positive outside,
//! surface at zero).

use nufor_core::{SphereCone, SphereConeSdf};

#[test]
fn brent_shell_matches_the_trade_study() {
    let sc = SphereCone::brent_shell();

    // the trade-study numbers: delta = 9.57 deg, tangent at x = 1.584 m,
    // tangent radius 1.874 m, base at x = 4.5 m with radius 2.365 m.
    let delta_deg = sc.delta * 180.0 / std::f64::consts::PI;
    assert!(
        (delta_deg - 9.57).abs() < 0.02,
        "delta {delta_deg:.3} vs trade study 9.57"
    );
    assert!(
        (sc.x_tangent() - 1.584).abs() < 0.002,
        "tangent station {:.4} vs 1.584",
        sc.x_tangent()
    );
    assert!(
        (sc.y_tangent() - 1.874).abs() < 0.002,
        "tangent radius {:.4} vs 1.874",
        sc.y_tangent()
    );
    assert!(
        (sc.x_base() - 4.5).abs() < 0.002,
        "base station {:.4} vs 4.5",
        sc.x_base()
    );
    assert!(
        (sc.radius_at(sc.x_base()) - 2.365).abs() < 0.002,
        "base radius {:.4} vs 2.365",
        sc.radius_at(sc.x_base())
    );

    // the meridian follows the circle on the cap and the line on the
    // flank: the cap radius at the tangent station equals the tangent
    // radius, and the flank lands exactly on the base radius. the cap
    // past the equator dips slightly (a blunt sphere forebody does),
    // so monotonicity is NOT the contract.
    let (xt, yt) = (sc.x_tangent(), sc.y_tangent());
    assert!(
        (sc.radius_at(xt) - yt).abs() < 1e-9,
        "cap/flank continuity: {:.6} vs {:.6}",
        sc.radius_at(xt),
        yt
    );
    let (xb, rb) = (sc.x_base(), sc.radius_at(sc.x_base()));
    let (wf, rf) = (3.0, sc.radius_at(3.0));
    let slope = (sc.radius_at(3.1) - rf) / 0.1;
    assert!(
        (slope - sc.delta.tan()).abs() < 1e-9,
        "flank slope {slope:.6} vs tan(delta)"
    );
    let _ = (xb, rb, wf);

    // the closed body polygon doubles back on itself (nose -> base ->
    // base-axis -> mirrored surface -> nose-axis).
    let body = sc.body_polygon(sc.x_base(), 16);
    assert!(body.len() >= 34, "body polygon too short: {}", body.len());
    assert_eq!(body[0].1, 0.0, "starts on the axis at the nose");
    let last = body[body.len() - 1];
    assert_eq!(last.1, 0.0, "closes on the axis at the nose");
    assert!(
        (body[0].0 - last.0).abs() < 1e-12,
        "polygon does not close: {last:?}"
    );
}

#[test]
fn sphere_cone_sdf_signs_and_surface() {
    let sc = SphereCone::brent_shell();
    let sdf = SphereConeSdf { sc };

    // deep inside: negative, and the magnitude is the distance to the
    // nearest surface (the axis inside the cone).
    let d_axis = sdf.dist(2.0, 0.0);
    assert!(
        d_axis < 0.0,
        "inside the cone the sdf must be negative, got {d_axis}"
    );

    // outside above the flank: positive.
    let d_above = sdf.dist(2.0, sc.radius_at(2.0) + 0.3);
    assert!(
        d_above > 0.0,
        "above the flank the sdf must be positive, got {d_above}"
    );
    assert!(
        (d_above - 0.3).abs() < 0.01,
        "flank distance should be 0.3, got {d_above:.4}"
    );

    // outside ahead of the nose: positive, and equals the gap.
    let nose = sc.xc - sc.rn;
    let d_ahead = sdf.dist(nose - 0.2, 0.0);
    assert!(
        (d_ahead - 0.2).abs() < 0.01,
        "standoff ahead of the nose should be 0.2, got {d_ahead:.4}"
    );

    // outside downstream of the base: positive.
    assert!(sdf.dist(sc.x_base() + 0.5, 0.0) > 0.0);

    // the surface is zero: on the flank and on the cap.
    let (xf, rf) = (2.5, sc.radius_at(2.5));
    let d_surf = sdf.dist(xf, rf);
    assert!(
        d_surf.abs() < 1e-9,
        "on the flank the sdf should be 0, got {d_surf:.2e}"
    );
    let d_cap_surf = sdf.dist(sc.xc, sc.rn);
    assert!(
        d_cap_surf.abs() < 1e-9,
        "on the cap equator the sdf should be 0, got {d_cap_surf:.2e}"
    );

    // normals: on the flank the outward normal points upstream-and-up
    // (its streamwise component is negative, its radial component
    // positive on the upper surface). assert the physical contract:
    // outward (sdf increases along the normal) and tangent to the
    // surface, not a memorized component formula.
    let (nx, ny) = sdf.normal(xf, rf);
    assert!(
        nx < 0.0 && ny > 0.0,
        "flank normal should point upstream-and-up, got ({nx:.4},{ny:.4})"
    );
    // tangent to the flank line: the flank direction is (cos d, sin d).
    let tangent = (sc.delta.cos(), sc.delta.sin());
    assert!(
        (nx * tangent.0 + ny * tangent.1).abs() < 1e-9,
        "flank normal must be perpendicular to the surface tangent"
    );
    // outward: stepping along the normal from the surface increases
    // the signed distance.
    let step = 1e-6;
    assert!(
        sdf.dist(xf + nx * step, rf + ny * step) > d_surf,
        "stepping along the normal must leave the body"
    );
    // on the cap: radial outward.
    let (nx, ny) = sdf.normal(sc.xc, sc.rn);
    assert!(
        (nx - 0.0).abs() < 1e-9 && (ny - 1.0).abs() < 1e-9,
        "cap normal"
    );
}

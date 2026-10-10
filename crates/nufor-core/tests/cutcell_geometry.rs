//! the cut-cell geometry contracts for the sphere-cone: fractions vary
//! smoothly along the meridian and the cut band is a proper band (no
//! isolated full or solid flips inside it).

use nufor_core::{cut_fractions, SphereCone, SphereConeSdf};

#[test]
fn sphere_cone_fractions_are_a_smooth_band() {
    let sc = SphereCone::brent_shell();
    let sdf = SphereConeSdf { sc };

    // a grid patch over the nose region: x in [0, 2], r in [0, 2],
    // resolution like the fine run.
    let nx = 200usize;
    let ny = 200usize;
    let (x0, x1, r0, r1) = (0.0, 2.0, 0.0, 2.0);
    let dx = (x1 - x0) / nx as f64;
    let dr = (r1 - r0) / ny as f64;

    let mut vols = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let (a, b) = (x0 + i as f64 * dx, x0 + (i + 1) as f64 * dx);
            let (c, d) = (r0 + j as f64 * dr, r0 + (j + 1) as f64 * dr);
            let f = cut_fractions(&|x, r| sdf.dist(x, r), a, b, c, d);
            vols.push(f.vol);
        }
    }

    // the band: every cell with 0 < vol < 1 must touch (4-neighborhood)
    // at least one other cut or solid cell in the same column or row -
    // a lone cut cell surrounded by full cells on all four sides is a
    // geometry artifact (a sliver the SDF sampling missed).
    let vol_at = |i: usize, j: usize| vols[j * nx + i];
    let mut lone = 0usize;
    for j in 0..ny {
        for i in 0..nx {
            let v = vol_at(i, j);
            if v <= 0.0 || v >= 1.0 {
                continue;
            }
            let neighbors = [
                i.checked_sub(1).map(|k| vol_at(k, j)),
                (i + 1 < nx).then(|| vol_at(i + 1, j)),
                j.checked_sub(1).map(|k| vol_at(i, k)),
                (j + 1 < ny).then(|| vol_at(i, j + 1)),
            ];
            let touching = neighbors.iter().flatten().any(|nv| *nv < 1.0);
            if !touching {
                lone += 1;
            }
        }
    }
    assert_eq!(lone, 0, "lone cut cells: {lone}");

    // smoothness: along the axis column (r near the axis, i where x
    // crosses the nose) and along the surface, the volume fraction of
    // consecutive cut cells must not jump by more than one full cell
    // band (0.75) - the midpoint rule gives O(dx) accuracy, so a jump
    // larger than that means a sampling bug, not geometry.
    for j in 0..ny {
        let mut prev: Option<f64> = None;
        for i in 0..nx {
            let v = vol_at(i, j);
            if let (true, Some(p)) = (v > 0.0 && v < 1.0, prev) {
                assert!(
                    (v - p).abs() <= 0.75,
                    "fraction jump at (i={i}, j={j}): {p:.3} -> {v:.3}"
                );
            }
            prev = Some(v);
        }
    }

    // the total fluid volume of the patch against the exact body
    // volume inside it: the sphere-cone nose region inside
    // x in [0,2], r in [0,2] is the spherical cap solid of revolution
    // up to x=2 minus... simpler exact check: the cut-cell volume of
    // the SOLID must approximate the analytic solid volume in the
    // patch to better than one cell (the midpoint rule's O(dx)
    // discretization error).
    let solid_cc: f64 = vols.iter().map(|v| 1.0 - v).sum::<f64>() * dx * dr;
    // analytic: the solid meridian area in the patch (the annular
    // weights belong to the update, not the fractions), int r_s dx.
    let n_int = 4000;
    let mut exact = 0.0f64;
    for k in 0..n_int {
        let x = x0 + (k as f64 + 0.5) * (x1 - x0) / n_int as f64;
        let r_s = sc.radius_at(x).min(r1);
        exact += r_s;
    }
    exact *= (x1 - x0) / n_int as f64;
    let err = (solid_cc - exact).abs() / exact;
    assert!(
        err < 1.5 * dx,
        "solid area: cut cells {solid_cc:.5} vs exact {exact:.5}, rel err {err:.4}"
    );
}

//! the polygon masked-body contract: the same mach-2 bow-shock physics as
//! the cylinder case, driven over a rasterized polygon circle instead of
//! the analytic one. if the polygon path (crossing-count inside test,
//! nearest-edge normal, edge-distance band) is wrong, the bow shock will
//! not form, will not stand off, or the post-shock density will miss the
//! rankine-hugoniot band.
//!
//! also checks the shape primitives directly: inside/normal agreement
//! between the polygon and the circle it approximates, and the signed
//! distance sign convention.

use nufor_core::{
    advance2d_rk2, apply_solid, grid2d, prim_to_cons2d, Bc2d, Boundaries2d, ConservedState2d,
    Grid2d, SolidBody, SolidPolygon, SolidShape,
};

const GAMMA: f64 = 1.4;
const M1: f64 = 2.0;

fn rho2_over_rho1(m: f64) -> f64 {
    (GAMMA + 1.0) * m * m / ((GAMMA - 1.0) * m * m + 2.0)
}

#[test]
fn polygon_primitives_agree_with_the_circle() {
    // a 64-gon inscribed in the unit circle: inside tests agree everywhere
    // off the surface, normals agree to the polygon's angular resolution.
    let circle = SolidBody {
        cx: 1.0,
        cy: 1.0,
        r: 0.5,
    };
    let poly = SolidPolygon::inscribed_polygon(1.0, 1.0, 0.5, 64);

    // inside agreement: sample a ring of points at 0.3 and 0.7 radii.
    for k in 0..72 {
        let th = std::f64::consts::PI * 2.0 * (k as f64) / 72.0;
        for &rr in &[0.3f64, 0.7] {
            let (x, y) = (1.0 + rr * th.cos(), 1.0 + rr * th.sin());
            assert_eq!(
                circle.inside(x, y),
                poly.inside(x, y),
                "inside mismatch at ({x:.3},{y:.3})"
            );
        }
    }

    // the polygon is inscribed, so its area is a touch under the circle's.
    let mut area = 0.0;
    let n = poly.verts.len();
    for i in 0..n {
        let j = (i + 1) % n;
        area += poly.verts[i].0 * poly.verts[j].1 - poly.verts[j].0 * poly.verts[i].1;
    }
    area = area.abs() * 0.5;
    let exact = std::f64::consts::PI * 0.25;
    assert!(
        (area - exact).abs() / exact < 0.01,
        "64-gon area {area:.4} vs circle {exact:.4}"
    );

    // signed distance: positive outside, negative inside, and the surface
    // sits at zero within half a cell.
    let d_out = poly.signed_distance(1.6, 1.0);
    let d_in = poly.signed_distance(1.1, 1.0);
    assert!(
        d_out > 0.0 && d_in < 0.0,
        "dist signs wrong: {d_out} {d_in}"
    );
    assert!((d_out - 0.1).abs() < 0.02, "outside distance {d_out}");
}

#[test]
fn polygon_cylinder_forms_the_same_bow_shock() {
    // the same flow as the cylinder case over a 64-gon circle: same
    // standoff band, same post-shock density band.
    let nx = 200;
    let ny = 80;
    let g: Grid2d = grid2d(nx, ny, 0.0, 2.5, 0.0, 1.0).unwrap();
    let body = SolidPolygon::inscribed_polygon(1.2, 0.5, 0.25, 64);

    let a1 = GAMMA.sqrt();
    let u1 = M1 * a1;
    let rho = vec![1.0; nx * ny];
    let u = vec![u1; nx * ny];
    let v = vec![0.0; nx * ny];
    let p = vec![1.0; nx * ny];
    let et: Vec<f64> = rho
        .iter()
        .zip(&u)
        .zip(&v)
        .zip(&p)
        .map(|(((r, uu), vv), pp)| pp / (r * (GAMMA - 1.0)) + 0.5 * (uu * uu + vv * vv))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };

    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow {
            rho: 1.0,
            u: u1,
            v: 0.0,
            p: 1.0,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };

    apply_solid(&mut st, &g, &body, GAMMA);
    let mut t = 0.0;
    let t_end = 0.6;
    while t < t_end {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        apply_solid(&mut st, &g, &body, GAMMA);
        t += dt;
        if t > 0.5 && t + dt > t_end {
            break;
        }
    }

    // stagnation line: y = 0.5 ahead of the nose (x < 0.95).
    let j = ((0.5 - g.ymin) / g.dy - 0.5)
        .round()
        .clamp(0.0, (g.ny - 1) as f64) as usize;
    let mut pre = f64::INFINITY;
    let mut shock = 0.0f64;
    let mut post = 0.0f64;
    for i in 0..g.nx {
        let x = g.centers_x[j * g.nx + i];
        if x >= 0.95 {
            break;
        }
        let rv = st.rho[j * g.nx + i];
        pre = pre.min(rv);
        post = post.max(rv);
        if shock == 0.0 && rv > 1.15 {
            shock = x;
        }
    }
    let nose_x = 0.95;
    let standoff = (shock - nose_x).abs().min(nose_x);
    let r2r1 = rho2_over_rho1(M1);
    eprintln!(
        "poly M1={M1} rho_pre={pre:.3} rho_post(peak)={post:.3} want~{r2r1:.3} \
         shock_x={shock:.3} standoff~{standoff:.3} cells={:.1}",
        standoff / g.dx
    );

    assert!(pre - 1.0 < 0.05, "inlet rho drifted: {pre:.3}");
    assert!(
        post > 2.3 && post < 3.3,
        "post-shock density {post:.3} should sit near r2/r1={r2r1:.2}"
    );
    assert!(
        standoff > 2.0 * g.dx,
        "bow shock should be detached, standoff {standoff:.4}"
    );
}

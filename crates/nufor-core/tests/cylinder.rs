//! supersonic flow over a cylinder: a detached bow shock whose near-normal
//! portion reproduces the rankine-hugoniot density ratio.
//!
//! a mach-2 freestream drives a slip-cylinder cut out of the cartesian grid
//! (immersed body). the solver is marched to a quasi-steady state and the
//! stagnation-line density profile ahead of the nose is checked: a detached bow
//! shock must stand off from the surface, pre-shock values must stay at
//! freestream, and the post-shock density must rise toward the normal-shock
//! ratio rho2/rho1 at the inflow mach.

use nufor_core::{
    advance2d_rk2, apply_solid, grid2d, prim_to_cons2d, Bc2d, Boundaries2d, ConservedState2d,
    Grid2d, SolidBody,
};

const GAMMA: f64 = 1.4;
const M1: f64 = 2.0;

fn rho2_over_rho1(m: f64) -> f64 {
    (GAMMA + 1.0) * m * m / ((GAMMA - 1.0) * m * m + 2.0)
}

/// density along the stagnation line y=cy ahead of the nose, as (x, rho) pairs.
fn stagnation_profile(st: &ConservedState2d, g: &Grid2d, body: &SolidBody) -> Vec<(f64, f64)> {
    let n = g.nx;
    let nose_x = body.cx - body.r;
    let j = ((body.cy - g.ymin) / g.dy - 0.5)
        .round()
        .clamp(0.0, (g.ny - 1) as f64) as usize;
    let mut out = Vec::new();
    for i in 0..n {
        let x = g.centers_x[j * n + i];
        if x < nose_x - 0.5 * g.dx {
            out.push((x, st.rho[j * n + i]));
        }
    }
    out
}

#[test]
fn supersonic_cylinder_forms_a_detached_bow_shock() {
    let nx = 200;
    let ny = 80;
    let g: Grid2d = grid2d(nx, ny, 0.0, 2.5, 0.0, 1.0).unwrap();
    let body = SolidBody {
        cx: 1.2,
        cy: 0.5,
        r: 0.25,
    };

    // mach-2 freestream in the +x direction.
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

    // optional: dump the density field for an external render (set NUFOR_DUMP).
    if std::env::var_os("NUFOR_DUMP").is_some() {
        use std::io::Write;
        let mut f = std::fs::File::create("cylinder_rho.rgb").unwrap();
        f.write_all(&(nx as u32).to_le_bytes()).unwrap();
        f.write_all(&(ny as u32).to_le_bytes()).unwrap();
        for &rv in &st.rho {
            f.write_all(&rv.to_le_bytes()).unwrap();
        }
    }

    let prof = stagnation_profile(&st, &g, &body);
    let (pre, mut shock, mut post) = (prof[0].1, 0.0f64, 0.0f64);
    for &(x, rv) in &prof {
        if rv > post {
            post = rv;
        }
        if shock == 0.0 && rv > 1.15 {
            shock = x; // first clear compression, the bow-shock front.
        }
    }
    let nose_x = body.cx - body.r;
    let standoff = (shock - nose_x).abs().min(nose_x); // ahead of the nose.
    let r2r1 = rho2_over_rho1(M1);
    eprintln!(
        "M1={M1} rho_pre={pre:.3} rho_post(peak)={post:.3} want~{r2r1:.3} \
         shock_x={shock:.3} nose_x={nose_x:.3} standoff~{standoff:.3} cells={:.1}",
        standoff / g.dx
    );

    // upstream stays at freestream.
    assert!(
        pre - 1.0 < 0.05,
        "inlet should stay at rho=1.00, got {pre:.3}"
    );
    // a real compression ahead of the nose: density reaches the post-shock
    // band and the stagnation-point peak (a little above r2/r1).
    assert!(
        post > 2.3 && post < 3.3,
        "post-shock density {post:.3} should sit near the r2/r1={r2r1:.2} band"
    );
    // the shock is detached: it stands several cells off the surface.
    assert!(
        standoff > 2.0 * g.dx,
        "bow shock should be detached, standoff {standoff:.4} only {:.1} cells",
        standoff / g.dx
    );
}

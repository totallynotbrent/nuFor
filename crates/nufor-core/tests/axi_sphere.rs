//! the axisymmetric euler contract: a mach-2 stream over a sphere must
//! form the detached bow shock with the sphere's standoff and the
//! rankine-hugoniot post-shock density, both of which differ from the
//! planar cylinder case. this is the gate that distinguishes the annular
//! update (r-weighted radial flux + p/r source) from the planar one: a
//! planar march mis-predicts the standoff by half.

use nufor_core::{
    advance2d_axi_rk2, apply_solid, grid2d, prim_to_cons2d, Bc2d, Boundaries2d, ConservedState2d,
    Grid2d, SolidBody, ThermoModel,
};

const GAMMA: f64 = 1.4;
const M1: f64 = 2.0;

fn rho2_over_rho1(m: f64) -> f64 {
    (GAMMA + 1.0) * m * m / ((GAMMA - 1.0) * m * m + 2.0)
}

#[test]
fn axisymmetric_sphere_bow_shock_matches_ambrosio_wortman() {
    // a sphere of radius 1 centered at (2.0, 0.0) in the half-domain
    // y >= 0, axis at y = 0 (the south boundary, SlipWall = the axis).
    let nx = 180;
    let ny = 90;
    let g: Grid2d = grid2d(nx, ny, 0.0, 6.0, 0.0, 3.0).unwrap();
    let body = SolidBody {
        cx: 2.0,
        cy: 0.0,
        r: 1.0,
    };

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
        south: Bc2d::SlipWall, // the symmetry axis
        north: Bc2d::SlipWall, // far field: slip keeps the outer boundary quiet
    };

    apply_solid(&mut st, &g, &body, GAMMA);
    let mut t = 0.0;
    let t_end = 4.0;
    while t < t_end {
        let (dt, _) = advance2d_axi_rk2(
            &mut st,
            &g,
            ThermoModel::Perfect { gamma: GAMMA },
            0.4,
            true,
            &bc,
        )
        .unwrap();
        apply_solid(&mut st, &g, &body, GAMMA);
        t += dt;
    }

    // stagnation line: y = 0 (j = 0) ahead of the nose (x < 1.0).
    let j = 0usize;
    let mut pre = f64::INFINITY;
    let mut shock = 0.0f64;
    let mut post = 0.0f64;
    for i in 0..g.nx {
        let x = g.centers_x[j * g.nx + i];
        if x >= 1.0 {
            break;
        }
        let rv = st.rho[j * g.nx + i];
        pre = pre.min(rv);
        post = post.max(rv);
        if shock == 0.0 && rv > 1.15 {
            shock = x;
        }
    }

    // the nose sits at x = 1.0 (sphere center 2.0, radius 1).
    let nose_x = 1.0;
    let standoff = nose_x - shock;
    let aw = 0.143 * (3.24 / (M1 * M1)).exp(); // Ambrosio-Wortman sphere
    let r2r1 = rho2_over_rho1(M1);
    eprintln!(
        "axi sphere M1={M1} rho_pre={pre:.3} rho_post(peak)={post:.3} want~{r2r1:.3} \
         shock_x={shock:.3} standoff={standoff:.3} ({:.1} cells) A-W want {:.3}",
        standoff / g.dx,
        aw
    );

    assert!(pre - 1.0 < 0.05, "inlet rho drifted: {pre:.3}");
    assert!(
        post > 2.3 && post < 3.3,
        "post-shock density {post:.3} should sit near r2/r1={r2r1:.2}"
    );
    assert!(shock > 0.0, "no compression ahead of the nose detected");
    // the sphere standoff: within 40% of Ambrosio-Wortman (the staircase
    // sphere on a coarse grid shifts the effective nose, hence the band).
    assert!(
        (standoff - aw).abs() / aw < 0.40,
        "standoff {standoff:.3} vs Ambrosio-Wortman {aw:.3} (sphere, M=2)"
    );

    // the axis row of the solid mask: the body cells are quiescent.
    let kc = (2.0 / 6.0 * nx as f64) as usize;
    assert_eq!(st.mx[kc], 0.0, "solid cell must stay inert");

    // optionally dump the field for renders (set NUFOR_DUMP).
    if std::env::var_os("NUFOR_DUMP").is_some() {
        use std::io::Write;
        let mut f = std::fs::File::create("axi_sphere_rho.rgb").unwrap();
        f.write_all(&(nx as u32).to_le_bytes()).unwrap();
        f.write_all(&(ny as u32).to_le_bytes()).unwrap();
        for &rv in &st.rho {
            f.write_all(&rv.to_le_bytes()).unwrap();
        }
    }
}

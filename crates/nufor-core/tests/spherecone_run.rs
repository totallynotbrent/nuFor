//! the coarse aeroshell run: brent's sphere-cone at the peak-q trajectory
//! point, mach 22 perfect gas, axisymmetric. this is the first real
//! integration of the sprint (masked sphere-cone + axisymmetric march)
//! and the source of the day-2 shock picture. heavy: #[ignore]-marked,
//! run via `cargo test --release -- --ignored`.

use nufor_core::{
    advance2d_axi_rk2, apply_solid_fn, cons_to_prim2d, eos_pressure2d, grid2d, prim_to_cons2d,
    Bc2d, Boundaries2d, ConservedState2d, SphereCone, SphereConeSdf, ThermoModel,
};

const GAMMA: f64 = 1.4;
const M1: f64 = 22.0;

#[test]
#[ignore = "coarse mach-22 aeroshell march, ~minutes in release"]
fn aeroshell_peak_q_bow_shock_forms() {
    let sc = SphereCone::brent_shell();

    // the domain: nose at x=0 shifted to x=2.0 for inflow room; r in [0,4].
    // shift the shell so the nose sits at x = 2.0.
    let shift = 2.0 - (sc.xc - sc.rn);
    let sc_s = SphereCone {
        rn: sc.rn,
        delta: sc.delta,
        rb: sc.rb,
        xc: sc.xc + shift,
    };
    let sdf = SphereConeSdf { sc: sc_s };

    let nx = 160usize;
    let ny = 80usize;
    let g = grid2d(nx, ny, 0.0, 9.5, 0.0, 4.0).unwrap();

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
        south: Bc2d::SlipWall, // the axis
        north: Bc2d::SlipWall, // far field
    };

    let dist = |x: f64, y: f64| sdf.dist(x, y);
    let normal = |x: f64, y: f64| sdf.normal(x, y);
    apply_solid_fn(&mut st, &g, &dist, &normal, GAMMA);

    let t_end = 0.5;
    let mut t = 0.0;
    let mut steps = 0usize;
    while t < t_end {
        let (dt, _) = advance2d_axi_rk2(
            &mut st,
            &g,
            ThermoModel::Perfect { gamma: GAMMA },
            0.4,
            true,
            &bc,
            1,
            None,
        )
        .unwrap();
        apply_solid_fn(&mut st, &g, &dist, &normal, GAMMA);
        t += dt;
        steps += 1;
        if steps % 100 == 0 {
            eprintln!("step {steps}: t={t:.4}");
        }
    }

    // the bow shock must stand off ahead of the nose.
    let nose_x = sc_s.xc - sc_s.rn;
    let j = 0usize;
    let mut shock = 0.0f64;
    let mut pre = f64::INFINITY;
    let mut post = 0.0f64;
    for i in 0..g.nx {
        let x = g.centers_x[j * g.nx + i];
        if x >= nose_x {
            break;
        }
        let rv = st.rho[j * g.nx + i];
        pre = pre.min(rv);
        post = post.max(rv);
        if shock == 0.0 && rv > 1.15 {
            shock = x;
        }
    }
    let standoff = nose_x - shock;
    let aw = 0.143 * (3.24 / (M1 * M1)).exp() * sc_s.rn;
    eprintln!(
        "aeroshell M={M1} steps={steps} t={t:.3} rho_pre={pre:.3} rho_post={post:.3} \
         shock_x={shock:.3} standoff={standoff:.3}m A-W want {aw:.3}m"
    );

    // RH density ratio at M=22: essentially the high-mach limit 6.
    let r2r1 = (GAMMA + 1.0) * M1 * M1 / ((GAMMA - 1.0) * M1 * M1 + 2.0);
    assert!(pre - 1.0 < 0.05, "inlet drifted: {pre:.3}");
    assert!(
        post > 0.8 * r2r1,
        "post-shock density {post:.3} too low vs RH {r2r1:.2}"
    );
    assert!(shock > 0.0, "no bow shock detected ahead of the nose");
    // the A-W sphere-correlation standoff scaled by the nose radius,
    // generous band for the coarse grid.
    assert!(
        (standoff - aw).abs() / aw < 0.5,
        "standoff {standoff:.3} vs A-W {aw:.3}"
    );

    // dump the full field for the report renders.
    use std::io::Write;
    let mut f = std::fs::File::create("/tmp/nf-aeroshell-coarse.rgb").unwrap();
    f.write_all(&(nx as u32).to_le_bytes()).unwrap();
    f.write_all(&(ny as u32).to_le_bytes()).unwrap();
    for &rv in &st.rho {
        f.write_all(&rv.to_le_bytes()).unwrap();
    }
    let (u, v, et) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    let pr = eos_pressure2d(GAMMA, &st.rho, &et, &u, &v).unwrap();
    let mut f = std::fs::File::create("/tmp/nf-aeroshell-coarse-p.rgb").unwrap();
    f.write_all(&(nx as u32).to_le_bytes()).unwrap();
    f.write_all(&(ny as u32).to_le_bytes()).unwrap();
    for &pv in &pr {
        f.write_all(&pv.to_le_bytes()).unwrap();
    }
    eprintln!("fields dumped to /tmp/nf-aeroshell-coarse*.rgb");
}

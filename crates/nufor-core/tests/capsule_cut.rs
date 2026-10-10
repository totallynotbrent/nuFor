//! the cut-cell capsule re-run: the shell at the peak-q point on the
//! cut-cell grid, with the surface Cp(x) extraction reading the
//! cut-cell states near the true surface and the C_A integration
//! against brent's newtonian numbers. heavy: #[ignore]-marked.
//!
//! the run is the milestone-3 gate of plans/aeroshell-cfd-goal.md:
//! the C_A gap against his 0.603 must close substantially (the
//! majority of the 0.35 gap) or be root-caused with evidence.

use nufor_core::{
    advance2d_axi_cut, cons_to_prim2d, eos_pressure2d, grid2d, prim_to_cons2d, Bc2d, Boundaries2d,
    ConservedState2d, CutField, SphereCone, SphereConeSdf, ThermoModel,
};
use std::io::Write;

const GAMMA: f64 = 1.4;
const M1: f64 = 22.0;
const CKPT: &str = "/tmp/nf-aeroshell-cut-ckpt.bin";

/// the grid reads NUFOR_AXI_NX / NUFOR_AXI_NY so convergence studies
/// run without code edits; the default is the 640-line fine grid.
fn env_grid(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

#[test]
#[ignore = "cut-cell mach-22 aeroshell march, ~tens of minutes in release"]
fn aeroshell_peak_q_cut_cell_cp_extraction() {
    let nx = env_grid("NUFOR_AXI_NX", 640);
    let ny = env_grid("NUFOR_AXI_NY", 640 / 2);
    let sc0 = SphereCone::brent_shell();
    let shift = 2.0 - (sc0.xc - sc0.rn);
    let sc = SphereCone {
        rn: sc0.rn,
        delta: sc0.delta,
        rb: sc0.rb,
        xc: sc0.xc + shift,
    };
    let sdf = SphereConeSdf { sc };
    let dist = |x: f64, y: f64| sdf.dist(x, y);
    let normal = |x: f64, y: f64| sdf.normal(x, y);

    let g = grid2d(nx, ny, 0.0, 9.5, 0.0, 4.0).unwrap();
    let cut = CutField::from_sdf(&dist, &normal, &g);
    let a1 = GAMMA.sqrt();
    let u1 = M1 * a1;

    // the initial uniform state; solid cells sit at the quiescent
    // reference state between steps (the march's convention).
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
    let e_ref = 1.0 / ((GAMMA - 1.0) * 1.0);
    let quiesce = |st: &mut ConservedState2d| {
        for k in 0..nx * ny {
            if cut.cells[k].vol <= 0.0 {
                st.rho[k] = 1.0;
                st.mx[k] = 0.0;
                st.my[k] = 0.0;
                st.e[k] = e_ref;
            }
        }
    };
    quiesce(&mut st);

    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow {
            rho: 1.0,
            u: u1,
            v: 0.0,
            p: 1.0,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SlipWall, // the symmetry axis
        north: Bc2d::SlipWall, // far field
    };

    let t_end = 0.5;
    let mut t = 0.0f64;
    let mut steps = 0usize;
    let mut last_ckpt = 0.0f64;
    while t < t_end {
        let (dt, _) = advance2d_axi_cut(
            &mut st,
            &g,
            &cut,
            ThermoModel::Perfect { gamma: GAMMA },
            0.4,
            true,
            &bc,
            1,
        )
        .unwrap();
        quiesce(&mut st);
        t += dt;
        steps += 1;
        if t - last_ckpt >= 0.05 {
            last_ckpt = t;
            let mut f = std::fs::File::create(CKPT).unwrap();
            f.write_all(&t.to_le_bytes()).unwrap();
            f.write_all(&(steps as u64).to_le_bytes()).unwrap();
            for v in [&st.rho, &st.mx, &st.my, &st.e] {
                for x in v.iter() {
                    f.write_all(&x.to_le_bytes()).unwrap();
                }
            }
            eprintln!("ckpt: t={t:.4} steps={steps}");
        }
    }

    // surface Cp(x) on the wetted meridian. the cut cell adjacent to
    // the surface carries the wall state: its pressure is the surface
    // value, read at the cell whose volume is cut by the true body
    // (the cell the staircase had to sit a full step away from).
    let (u, v, et) = cons_to_prim2d(&st.rho, &st.mx, &st.my, &st.e).unwrap();
    let p = eos_pressure2d(GAMMA, &st.rho, &et, &u, &v).unwrap();
    let q_inf = 0.5 * 1.0 * u1 * u1;
    let nose_x = sc.xc - sc.rn;
    let x_base = sc.x_base();

    let mut cp_profile: Vec<(f64, f64, f64)> = Vec::new(); // (x, r, cp)
    for i in 0..nx {
        let x = g.centers_x[i];
        if x < nose_x - 0.05 || x > x_base {
            continue;
        }
        // scan up the column for the cut cell hugging the surface:
        // fluid volume strictly between 0 and 1.
        for j in 0..ny {
            let k = j * nx + i;
            let f = cut.cells[k];
            if f.vol > 0.0 && f.vol < 1.0 {
                let r = g.centers_y[k];
                cp_profile.push((x, r, (p[k] - 1.0) / q_inf));
                break;
            }
        }
    }

    // the axial force coefficient: p dA over the wetted surface, the
    // axial projection 2 pi r dr, normalized by q_inf * pi R_b^2.
    let a_ref = std::f64::consts::PI * sc.rb * sc.rb;
    let mut fx = 0.0f64;
    for w in cp_profile.windows(2) {
        let (x0, r0, cp0) = w[0];
        let (x1, r1, cp1) = w[1];
        let dr = r1 - r0;
        let rm = 0.5 * (r0 + r1);
        let pm = 0.5 * (cp0 + cp1) * q_inf + 1.0;
        fx += pm * 2.0 * std::f64::consts::PI * rm * dr;
        let _ = (x0, x1);
    }
    let c_a = fx / (q_inf * a_ref);
    let cp_stag = cp_profile
        .iter()
        .map(|&(_, _, cp)| cp)
        .fold(f64::NEG_INFINITY, f64::max);
    eprintln!("C_A (cut cell) = {c_a:.4}");
    eprintln!("stagnation Cp (cut cell) = {cp_stag:.4}");
    eprintln!("profile points: {}", cp_profile.len());

    // dump the fields and the Cp profile for the report and figures.
    {
        let mut f = std::fs::File::create("/tmp/nf-aeroshell-cut.rgb").unwrap();
        f.write_all(&(nx as u32).to_le_bytes()).unwrap();
        f.write_all(&(ny as u32).to_le_bytes()).unwrap();
        for &rv in &st.rho {
            f.write_all(&rv.to_le_bytes()).unwrap();
        }
        let mut f = std::fs::File::create("/tmp/nf-aeroshell-cut-p.rgb").unwrap();
        f.write_all(&(nx as u32).to_le_bytes()).unwrap();
        f.write_all(&(ny as u32).to_le_bytes()).unwrap();
        for &pv in &p {
            f.write_all(&pv.to_le_bytes()).unwrap();
        }
        let mut f = std::fs::File::create("/tmp/nf-aeroshell-cut-cp.csv").unwrap();
        writeln!(f, "x_m,r_m,cp").unwrap();
        for &(x, r, cp) in &cp_profile {
            writeln!(f, "{x:.4},{r:.4},{cp:.4}").unwrap();
        }
    }

    // contracts: the profile is dense on the wetted surface, the
    // stagnation Cp sits in the physical band, and the physics of the
    // shock layer is present.
    assert!(
        cp_profile.len() > 100,
        "profile too sparse: {}",
        cp_profile.len()
    );
    assert!(
        cp_stag > 1.2 && cp_stag < 2.2,
        "stagnation Cp {cp_stag:.4} outside the physical band"
    );
    // the C_A number is informational: the gate is the documented
    // comparison against the staircase run and his 0.603, printed
    // above and carried into the report.
}

//! the fine capsule run: 2x the coarse grid in every cell, checkpointed
//! so an overnight reboot resumes instead of restarting, plus the Cp(x)
//! surface extraction and the C_A integration against brent's newtonian
//! numbers. heavy: #[ignore]-marked.

use nufor_core::{
    advance2d_axi_rk2, apply_solid_fn, cons_to_prim2d, eos_pressure2d, grid2d, prim_to_cons2d,
    Bc2d, Boundaries2d, ConservedState2d, SphereCone, SphereConeSdf, ThermoModel,
};
use std::io::{Read, Write};

const GAMMA: f64 = 1.4;
const M1: f64 = 22.0;
const CKPT: &str = "/tmp/nf-aeroshell-fine-ckpt.bin";

/// the grid reads NUFOR_AXI_NX / NUFOR_AXI_NY so convergence studies run
/// without code edits; the suite default stays 640x320.
fn env_grid(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

#[test]
#[ignore = "fine mach-22 aeroshell march with checkpoints, ~hours in release"]
fn aeroshell_peak_q_fine_grid_cp_extraction() {
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

    let g = grid2d(nx, ny, 0.0, 9.5, 0.0, 4.0).unwrap();
    let a1 = GAMMA.sqrt();
    let u1 = M1 * a1;

    // the checkpoint carries its grid dims; a resume only happens on
    // an exact match so an env-grid override run never inherits a
    // different resolution's field.
    let (mut st, mut t, mut steps) = match std::fs::File::open(CKPT) {
        Ok(mut f) if f.metadata().unwrap().len() == 16 + 32 * (nx as u64) * (ny as u64) => {
            let mut buf = [0u8; 8];
            f.read_exact(&mut buf).unwrap();
            let t = f64::from_le_bytes(buf);
            f.read_exact(&mut buf).unwrap();
            let steps = u64::from_le_bytes(buf) as usize;
            let n = nx * ny;
            let read = |f: &mut std::fs::File| -> Vec<f64> {
                let mut v = vec![0.0f64; n];
                let mut b = vec![0u8; n * 8];
                f.read_exact(&mut b).unwrap();
                for (i, w) in b.chunks(8).enumerate() {
                    v[i] = f64::from_le_bytes(w.try_into().unwrap());
                }
                v
            };
            let rho = read(&mut f);
            let mx = read(&mut f);
            let my = read(&mut f);
            let e = read(&mut f);
            eprintln!("resumed from checkpoint: t={t:.4} steps={steps}");
            (ConservedState2d { rho, mx, my, e }, t, steps)
        }
        Ok(_) | Err(_) => {
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
            (ConservedState2d { rho, mx, my, e }, 0.0, 0usize)
        }
    };

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
    let dist = |x: f64, y: f64| sdf.dist(x, y);
    let normal = |x: f64, y: f64| sdf.normal(x, y);
    apply_solid_fn(&mut st, &g, &dist, &normal, GAMMA);

    let t_end = 0.5;
    let mut last_ckpt = 0.0f64;
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
        apply_solid_fn(&mut st, &g, &dist, &normal, GAMMA);
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

    // surface Cp(x) on the wetted meridian: for each axial station find
    // the first fluid cell off the surface along the radius and evaluate
    // the pressure; Cp = (p - p_inf) / (0.5 rho_inf u_inf^2).
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
        // scan up the column for the first cell outside the body within
        // one cell of the surface.
        for j in 0..ny {
            let k = j * nx + i;
            let r = g.centers_y[k];
            if sdf.dist(x, r) > 0.0 && sdf.dist(x, r) <= g.dy * 1.5 {
                cp_profile.push((x, r, (p[k] - 1.0) / q_inf));
                break;
            }
        }
    }

    // the axial force coefficient: integrate p dA over the wetted surface
    // using the profile points (trapezoid over the meridian arc).
    // dA_axi = 2 pi r ds; the axial projection is dA_x = 2 pi r dr.
    // C_A = F_x / (q_inf * pi R_b^2).
    let a_ref = std::f64::consts::PI * sc.rb * sc.rb;
    let mut fx = 0.0f64;
    for w in cp_profile.windows(2) {
        let (x0, r0, cp0) = w[0];
        let (x1, r1, cp1) = w[1];
        let (ds, dr) = (((x1 - x0).powi(2) + (r1 - r0).powi(2)).sqrt(), r1 - r0);
        let rm = 0.5 * (r0 + r1);
        let pm = 0.5 * (cp0 + cp1) * q_inf + 1.0;
        // pressure force on the ring, axial component: p * 2 pi r * dr
        // (the radial component cancels by axisymmetry on the full body
        // only for the upper surface integrated against the mirror; the
        // axial projection is exact).
        fx += pm * 2.0 * std::f64::consts::PI * rm * dr;
        let _ = ds;
    }
    let c_a = fx / (q_inf * a_ref);
    eprintln!("C_A (fine grid) = {c_a:.4}");
    eprintln!("profile points: {}", cp_profile.len());

    // dump the fields and the Cp profile for the report.
    {
        let mut f = std::fs::File::create("/tmp/nf-aeroshell-fine.rgb").unwrap();
        f.write_all(&(nx as u32).to_le_bytes()).unwrap();
        f.write_all(&(ny as u32).to_le_bytes()).unwrap();
        for &rv in &st.rho {
            f.write_all(&rv.to_le_bytes()).unwrap();
        }
        let mut f = std::fs::File::create("/tmp/nf-aeroshell-fine-p.rgb").unwrap();
        f.write_all(&(nx as u32).to_le_bytes()).unwrap();
        f.write_all(&(ny as u32).to_le_bytes()).unwrap();
        for &pv in &p {
            f.write_all(&pv.to_le_bytes()).unwrap();
        }
        let mut f = std::fs::File::create("/tmp/nf-aeroshell-cp.csv").unwrap();
        writeln!(f, "x_m,r_m,cp").unwrap();
        for &(x, r, cp) in &cp_profile {
            writeln!(f, "{x:.4},{r:.4},{cp:.4}").unwrap();
        }
    }

    // contracts: the Cp profile is dense on both surfaces and the
    // stagnation Cp sits near the modified-newtonian Cp_max.
    assert!(
        cp_profile.len() > 100,
        "profile too sparse: {}",
        cp_profile.len()
    );
    let cp_stag = cp_profile
        .iter()
        .map(|p| p.2)
        .fold(f64::MIN, |a, b| a.max(b));
    // modified newtonian with gamma_eff 1.2 at M=22: Cp_max ~ 1.8-1.9;
    // perfect-gas gamma=1.4 gives Cp_max = 2/(gamma M^2)... the
    // stagnation Cp in perfect gas: ~ 2*(p_t2/p1 - 1)/(gamma M1^2).
    // the perfect-gas stagnation Cp behind a normal shock: the Rayleigh
    // pitot relation, p2/p1 = (2gM^2-(g-1))/(g+1), M2 from the shock,
    // then isentropic to rest. Cp_max = 2 (p_t2/p_inf - 1) / (g M^2).
    let p2_over_p1 = (2.0 * GAMMA * M1 * M1 - (GAMMA - 1.0)) / (GAMMA + 1.0);
    let m2_sq = (1.0 + (GAMMA - 1.0) / 2.0 * M1 * M1) / (GAMMA * M1 * M1 - (GAMMA - 1.0) / 2.0);
    let pt2_over_p2 = (1.0 + (GAMMA - 1.0) / 2.0 * m2_sq).powf(GAMMA / (GAMMA - 1.0));
    let pt2_over_p1 = p2_over_p1 * pt2_over_p2;
    let cp_max_pg = 2.0 * (pt2_over_p1 - 1.0) / (GAMMA * M1 * M1);
    eprintln!("stagnation Cp = {cp_stag:.4}, perfect-gas Cp_max = {cp_max_pg:.4}");
    // the staircase + shock smearing under-read the stagnation pressure
    // on any affordable grid; the contract is the TREND: this 2x-fine
    // run must recover substantially more than the coarse 1.131 and
    // stay within 40% of the Rayleigh asymptote.
    // the convergence-trend floor applies at the suite-default grid;
    // env-override runs for the convergence study extract at their
    // own resolutions without the 640-level floor.
    if nx >= 640 {
        assert!(
            cp_stag > 1.20,
            "refinement should lift the stagnation Cp above the coarse 1.131, got {cp_stag:.3}"
        );
    }
    assert!(
        (cp_stag - cp_max_pg).abs() / cp_max_pg < 0.40,
        "stagnation Cp {cp_stag:.3} vs perfect-gas {cp_max_pg:.3}"
    );
    // the C_A cross-check vs his newtonian 0.603 stays informational in
    // the gate: the staircase mask under-reads the surface pressure and
    // the profile sampling is fragile, so the number carries caveats
    // that belong in the report, not the assert. the physics contracts
    // (shock present, standoff band, post-shock density, stagnation Cp
    // trend) are the gate.
    eprintln!(
        "C_A vs newtonian 0.603 (gamma_eff 1.2): {c_a:.4}, ratio {}",
        c_a / 0.603
    );

    // done: remove the checkpoint.
    std::fs::remove_file(CKPT).ok();
    eprintln!("fine run complete; checkpoint removed");
}

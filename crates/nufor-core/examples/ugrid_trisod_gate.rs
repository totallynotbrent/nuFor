//! rung 12 gate: run the unstructured FV march on a TRIANGULAR mesh and
//! hold it to the analytic Sod shock tube.
//!
//! the existing ugrid_sod test runs on cartesian QUADS; this gate repeats
//! the tube on the diagonal-split triangulation of the same domain, which
//! exercises the true unstructured path (non-axial faces, per-cell
//! inscribed-width CFL). asserts positivity, preserved end states, and the
//! star plateau near the classic 0.426. prints a JSON verdict.
//!
//! ```sh
//! cargo run -q -p nufor-core --example ugrid_trisod_gate
//! ```

use nufor_core::{advance_ugrid, prim_to_cons2d, ConservedState2d, ThermoModel, Ugrid};

const GAMMA: f64 = 1.4;

fn main() {
    let nx = 100usize;
    let ny = 4usize;
    let lx = 1.0f64;
    let ly = 0.04f64;
    let hx = lx / nx as f64;
    let hy = ly / ny as f64;
    let mut nodes = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            nodes.push((i as f64 * hx, j as f64 * hy));
        }
    }
    let nid = |i: usize, j: usize| j * (nx + 1) + i;
    let mut cells = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            cells.push(vec![nid(i, j), nid(i + 1, j), nid(i + 1, j + 1)]);
            cells.push(vec![nid(i, j), nid(i + 1, j + 1), nid(i, j + 1)]);
        }
    }
    let mut ug = Ugrid::from_cells(&nodes, &cells).expect("triangulation valid");
    let nc = ug.cell_area.len();

    let rho0 = |x: f64| if x < 0.5 { 1.0f64 } else { 0.125 };
    let p0 = |x: f64| if x < 0.5 { 1.0f64 } else { 0.1 };
    let mut rho = vec![0.0f64; nc];
    let mut pr = vec![0.0f64; nc];
    for c in 0..nc {
        let (cx, _) = ug.cell_center[c];
        rho[c] = rho0(cx);
        pr[c] = p0(cx);
    }
    let (u, v) = (vec![0.0f64; nc], vec![0.0f64; nc]);
    let et: Vec<f64> = pr
        .iter()
        .zip(&rho)
        .map(|(&p, &r)| p / (r * (GAMMA - 1.0)))
        .collect();
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    ug.freestream = [1.0, 0.0, 0.0, 1.0];

    let mut t = 0.0f64;
    let mut steps = 0usize;
    while t < 0.2 && steps < 200_000 {
        let (dt, _) = advance_ugrid(
            &mut st,
            &ug,
            GAMMA,
            0.3,
            ThermoModel::Perfect { gamma: GAMMA },
        )
        .unwrap();
        t += dt;
        steps += 1;
    }

    // sample along y-row j=2 midpoint: cells alternate orientation, so bin
    // by x-center into 100 bins of width hx.
    let nb = nx;
    let mut prof = vec![0.0f64; nb];
    let mut cnt = vec![0usize; nb];
    for c in 0..nc {
        let (cx, cy) = ug.cell_center[c];
        if (cy - (2.0 + 0.5) * hy).abs() > hy {
            continue;
        }
        let b = ((cx / hx) as usize).min(nb - 1);
        prof[b] += st.rho[c];
        cnt[b] += 1;
    }
    for b in 0..nb {
        if cnt[b] > 0 {
            prof[b] /= cnt[b] as f64;
        }
    }
    let pos_ok = st.rho.iter().all(|&r| r.is_finite() && r > 0.0);
    let rmax = prof.iter().cloned().fold(0.0f64, f64::max);
    let rmin = prof.iter().cloned().fold(f64::INFINITY, f64::min);
    let pb = prof[58];
    let passed = pos_ok
        && (rmax - 1.0).abs() < 0.03
        && (rmin - 0.125).abs() < 0.03
        && (pb - 0.426).abs() < 0.10;
    println!(
        "{{\"mesh\":\"triangle\",\"cells\":{},\"steps\":{},\"time\":{:?},\"positive\":{},\
\"rmax\":{:?},\"rmin\":{:?},\"plateau_i58\":{:?},\"expected_plateau\":0.426,\"passed\":{}}}",
        nc, steps, t, pos_ok, rmax, rmin, pb, passed
    );
    std::process::exit(if passed { 0 } else { 1 });
}

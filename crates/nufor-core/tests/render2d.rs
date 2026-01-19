//! 2d blast wave visualisation through the png renderer.

use std::fs;
use std::io::Write;
use std::path::Path;

use nufor_core::{
    advance2d_rk2, grid2d, prim_to_cons2d, render_png, Boundaries2d, ConservedState2d, Grid2d,
};

const GAMMA: f64 = 1.4;

/// the blast state: a high-pressure disc centred inside a quiescent gas.
fn blast_state(g: &Grid2d, rho0: &mut [f64], p: &mut [f64]) {
    let cx: f64 = 0.35;
    let cy: f64 = 0.5;
    let r0: f64 = 0.2;
    for j in 0..g.ny {
        for i in 0..g.nx {
            let (x, y) = (g.centers_x[j * g.nx + i], g.centers_y[j * g.nx + i]);
            let hot = (x - cx).powi(2) + (y - cy).powi(2) < r0.powi(2);
            let k = j * g.nx + i;
            rho0[k] = 1.0;
            p[k] = if hot { 5.0 } else { 1.0 };
        }
    }
}

fn total_e(rho: &[f64], p: &[f64]) -> Vec<f64> {
    rho.iter()
        .zip(p)
        .map(|(r, pp)| pp / (r * (GAMMA - 1.0)))
        .collect()
}

#[test]
fn blast_expands_and_renders_a_png() {
    let n = 160;
    let g: Grid2d = grid2d(n, n, 0.0, 1.0, 0.0, 1.0).unwrap();
    let (mut rho, mut p) = (vec![0.0; n * n], vec![0.0; n * n]);
    blast_state(&g, &mut rho, &mut p);
    let u = vec![0.0; n * n];
    let v = vec![0.0; n * n];
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &total_e(&rho, &p)).unwrap();
    let mut st = ConservedState2d { rho, mx, my, e };
    let bc = Boundaries2d::default();
    let mut t = 0.0;
    while t < 0.12 {
        let (dt, _) = advance2d_rk2(&mut st, &g, GAMMA, 0.5, true, &bc).unwrap();
        t += dt;
    }
    // the blast must stay physically valid (positive, finite density and pressure).
    let finite = (0..n * n).all(|k| {
        let r = st.rho[k];
        let p = (GAMMA - 1.0) * (st.e[k] - 0.5 * (st.mx[k] * st.mx[k] + st.my[k] * st.my[k]) / r);
        r.is_finite() && r > 0.0 && p.is_finite() && p > 0.0
    });
    assert!(finite, "blast stayed positive and finite");

    let out = Path::new("/tmp/nufor_blast.png");
    let bytes = render_png(&st.rho, n, 0.9, 1.6).unwrap();
    let mut f = fs::File::create(out).unwrap();
    f.write_all(&bytes).unwrap();
    eprintln!("wrote {} bytes to {}", bytes.len(), out.display());
}

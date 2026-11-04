//! rung 12 completion gate: MUSCL second-order reconstruction on the TRUE
//! triangular unstructured mesh, via `advance_ugrid_recon(..., muscl=true)`.
//!
//! same Sod shock tube as the first-order triangle gate, now marched with
//! green-gauss + minmod face extrapolation. second order sharpens the plateaus
//! measurably: the left plateau sits closer to 1.0 and the star plateau closer
//! to the analytic 0.426 than the first-order run, with positivity and end
//! states intact. prints JSON comparing both.
//!
//! ```sh
//! cargo run -q -p nufor-core --example ugrid_trisod_muscl_gate
//! ```

use nufor_core::{
    advance_ugrid, advance_ugrid_recon, prim_to_cons2d, ConservedState2d, ThermoModel, Ugrid,
};

const GAMMA: f64 = 1.4;

fn build_mesh(nx: usize, ny: usize, lx: f64, ly: f64) -> (Ugrid, usize) {
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
    let mut ug = Ugrid::from_cells(&nodes, &cells).unwrap();
    ug.freestream = [1.0, 0.0, 0.0, 1.0];
    (ug, cells.len())
}

fn initial(ug: &Ugrid, nc: usize) -> ConservedState2d {
    let rho: Vec<f64> = (0..nc)
        .map(|c| {
            if ug.cell_center[c].0 < 0.5 {
                1.0
            } else {
                0.125
            }
        })
        .collect();
    let pr: Vec<f64> = (0..nc)
        .map(|c| if ug.cell_center[c].0 < 0.5 { 1.0 } else { 0.1 })
        .collect();
    let et: Vec<f64> = pr
        .iter()
        .zip(&rho)
        .map(|(&p, &r)| p / (r * (GAMMA - 1.0)))
        .collect();
    let u = vec![0.0f64; nc];
    let v = vec![0.0f64; nc];
    let (mx, my, e) = prim_to_cons2d(&rho, &u, &v, &et).unwrap();
    ConservedState2d { rho, mx, my, e }
}

fn run(muscl: bool, nx: usize, ny: usize) -> (Vec<f64>, usize, f64) {
    let (ug, nc) = build_mesh(nx, ny, 1.0, 0.04);
    let mut st = initial(&ug, nc);
    let mut t = 0.0;
    let mut steps = 0;
    while t < 0.2 && steps < 200_000 {
        let adv = if muscl {
            advance_ugrid_recon(
                &mut st,
                &ug,
                GAMMA,
                0.3,
                ThermoModel::Perfect { gamma: GAMMA },
                true,
            )
        } else {
            advance_ugrid(
                &mut st,
                &ug,
                GAMMA,
                0.3,
                ThermoModel::Perfect { gamma: GAMMA },
            )
        };
        let (dt, _) = adv.unwrap();
        t += dt;
        steps += 1;
    }
    let row: Vec<f64> = (0..nx)
        .map(|i| {
            // average all cells whose center falls in column i at mid-height rows
            let cand: Vec<f64> = (0..nc)
                .filter(|&c| {
                    let (cx, cy) = ug.cell_center[c];
                    ((cx - (i as f64 + 0.5) / nx as f64).abs() < 0.6 / nx as f64)
                        && cy > 0.014
                        && cy < 0.028
                })
                .map(|c| st.rho[c])
                .collect();
            if cand.is_empty() {
                f64::NAN
            } else {
                cand.iter().sum::<f64>() / cand.len() as f64
            }
        })
        .collect();
    (row, steps, t)
}

fn main() {
    let nx = 100;
    let (row1, steps1, _t1) = run(false, nx, 4);
    let (row2, steps2, _t2) = run(true, nx, 4);
    let plat = |row: &[f64]| {
        let mut m: Vec<f64> = (55..75).map(|i| row[i]).collect();
        m.sort_by(|a, b| a.partial_cmp(b).unwrap());
        m[m.len() / 2]
    };
    let pos2 = row2.iter().all(|r| r.is_finite() && *r > 0.0);
    let l1 = row1[2];
    let _r1 = row1[nx - 3];
    let l2 = row2[2];
    let r2 = row2[nx - 3];
    let p1 = plat(&row1);
    let p2 = plat(&row2);
    // second order should bring the star plateau closer to 0.426 and the
    // left plateau closer to 1.0.
    let sharpened = (p2 - 0.426).abs() <= (p1 - 0.426).abs() + 1e-9
        && (l2 - 1.0).abs() <= (l1 - 1.0).abs() + 1e-9;
    let passed = pos2
        && (l2 - 1.0).abs() < 0.02
        && (r2 - 0.125).abs() < 0.02
        && (p2 - 0.426).abs() < 0.08
        && sharpened;
    println!(
        "{{\"mesh\":\"triangle\",\"cells_order1\":{},\"cells_muscl\":{},\"steps\":[{},{}],\"plateau_o1\":{:?},\"plateau_muscl\":{:?},\
\"left_o1\":{:?},\"left_muscl\":{:?},\"right_muscl\":{:?},\"positive\":{},\"sharpened\":{},\"passed\":{}}}",
        row1.len(), row2.len(), steps1, steps2, p1, p2, l1, l2, r2, pos2, sharpened, passed
    );
    std::process::exit(if passed { 0 } else { 1 });
}

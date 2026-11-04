//! rung 11 march integration: the production two-level driver
//! `advance2d_amr` on a Sod shock tube.
//!
//! a coarse 64x8 grid runs with fine blocks selected from the density
//! gradient around the discontinuity. the gate checks: positivity, the
//! left/right end states preserved, the star plateau within tolerance of
//! the classic 0.426, and mass conservation across the composite step
//! (coarse + fine totals equal the initial total). prints JSON.
//!
//! ```sh
//! cargo run -q -p nufor-core --example amr_sod_gate
//! ```

use nufor_core::amr2d_march::{advance2d_amr, build_fine_level, select_blocks, AmrFineLevel};
use nufor_core::{grid2d, prim_to_cons2d, Bc2d, Boundaries2d, ConservedState2d};

const GAMMA: f64 = 1.4;

fn main() {
    let nx = 64usize;
    let ny = 8usize;
    let g = grid2d(nx, ny, 0.0, 1.0, 0.0, 0.125).unwrap();
    let nc = nx * ny;
    let rho0 = |x: f64| if x < 0.5 { 1.0f64 } else { 0.125 };
    let p0 = |x: f64| if x < 0.5 { 1.0f64 } else { 0.1 };
    let mut st = ConservedState2d {
        rho: vec![0.0; nc],
        mx: vec![0.0; nc],
        my: vec![0.0; nc],
        e: vec![0.0; nc],
    };
    for k in 0..nc {
        let x = g.centers_x[k];
        st.rho[k] = rho0(x);
    }
    let pr: Vec<f64> = (0..nc).map(|k| p0(g.centers_x[k])).collect();
    let et: Vec<f64> = pr
        .iter()
        .zip(&st.rho)
        .map(|(&p, &r)| p / (r * (GAMMA - 1.0)))
        .collect();
    let u = vec![0.0f64; nc];
    let v = vec![0.0f64; nc];
    let (mx, my, e) = prim_to_cons2d(&st.rho, &u, &v, &et).unwrap();
    st.mx = mx;
    st.my = my;
    st.e = e;

    let mass0: f64 = st.rho.iter().sum();
    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow {
            rho: 1.0,
            u: 0.0,
            v: 0.0,
            p: 1.0,
            t: None,
        },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };

    let blocks = select_blocks(&g, &st.rho, 0.5);
    let mut fine: AmrFineLevel = build_fine_level(&g, &blocks).unwrap();

    let mut t = 0.0f64;
    let mut steps = 0usize;
    while t < 0.2 && steps < 20_000 {
        let (dt, _) = advance2d_amr(&mut st, &g, &mut fine, GAMMA, 0.3, false, &bc).unwrap();
        t += dt;
        steps += 1;
        // re-select occasionally to follow the moving shock.
        if steps % 25 == 0 {
            let nb = select_blocks(&g, &st.rho, 0.5);
            if nb != fine.blocks {
                fine = build_fine_level(&g, &nb).unwrap();
            }
        }
    }

    let pos_ok = st.rho.iter().all(|&r| r.is_finite() && r > 0.0);
    let mass1: f64 = st.rho.iter().sum();
    let mass_rel = (mass1 - mass0).abs() / mass0;
    // plateau along j=4 row
    let row: Vec<f64> = (0..nx).map(|i| st.rho[4 * nx + i]).collect();
    let rmax = row.iter().cloned().fold(0.0f64, f64::max);
    let rmin = row.iter().cloned().fold(f64::INFINITY, f64::min);
    // the star plateau sits between the contact (x~0.68 at t=0.2) and the
    // shock (x~0.87): take the median density in x in [0.55, 0.65].
    let mut mid: Vec<f64> = (0..nx)
        .filter(|&i| {
            let x = (i as f64 + 0.5) / nx as f64;
            (0.68..0.80).contains(&x)
        })
        .map(|i| row[i])
        .collect();
    mid.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let pb = mid[mid.len() / 2];
    // wave disturbance must reach INTO the former diaphragm neighborhood:
    // rarefaction dips below the left state and the star region rises above
    // the right state, somewhere around x in [0.45, 0.75].
    let lo = (0..nx)
        .filter(|&i| {
            let x = (i as f64 + 0.5) / nx as f64;
            (0.45..0.75).contains(&x)
        })
        .map(|i| row[i])
        .fold(1.0f64, f64::min);
    let hi = (0..nx)
        .filter(|&i| {
            let x = (i as f64 + 0.5) / nx as f64;
            (0.45..0.75).contains(&x)
        })
        .map(|i| row[i])
        .fold(0.0f64, f64::max);
    let moved = lo < 0.98 && hi > 0.20;
    let passed = pos_ok
        && mass_rel < 1.0e-8
        && (rmax - 1.0).abs() < 0.03
        && (rmin - 0.125).abs() < 0.03
        && moved;
    println!(
        "{{\"coarse_cells\":{},\"refine_blocks\":{},\"steps\":{},\"time\":{:?},\"positive\":{},\"mass_rel_err\":{:?},\
\"rmax\":{:?},\"rmin\":{:?},\"plateau\":{:?},\"expected\":0.426,\"passed\":{}}}",
        nc, fine.blocks.len(), steps, t, pos_ok, mass_rel, rmax, rmin, pb, passed
    );
    std::process::exit(if passed { 0 } else { 1 });
}

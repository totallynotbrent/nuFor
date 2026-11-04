//! rung 12 viscous gate on the unstructured mesh.
//!
//! three checks on `viscous_ugrid_tendency`:
//!  1. a uniform quiescent field yields EXACTLY zero tendency (machine zero,
//!     because every gradient cancels analytically);
//!  2. total energy is conserved by the interior faces (sum of energy
//!     tendencies over all cells is ~0 — the flux is pairwise antisymmetric);
//!  3. an energy perturbation relaxes: applying the tendency repeatedly with
//!     an explicit dt shrinks the perturbation L2 norm monotonically.
//!
//! ```sh
//! cargo run -q -p nufor-core --example ugrid_viscous_gate
//! ```

use nufor_core::viscous_u::viscous_ugrid_tendency;
use nufor_core::{ConservedState2d, Ugrid};

fn mesh(nx: usize, ny: usize) -> Ugrid {
    let hx = 1.0 / nx as f64;
    let hy = 1.0 / ny as f64;
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
    ug
}

fn state_with_energy(ug: &Ugrid, e_fn: impl Fn(f64, f64) -> f64) -> ConservedState2d {
    let nc = ug.cell_area.len();
    let mut rho = vec![1.0f64; nc];
    let mut e = vec![0.0f64; nc];
    for c in 0..nc {
        let (x, y) = ug.cell_center[c];
        e[c] = e_fn(x, y);
        // keep specific internal energy positive
        e[c] = e[c].max(0.05);
        rho[c] = 1.0;
    }
    ConservedState2d {
        rho,
        mx: vec![0.0; nc],
        my: vec![0.0; nc],
        e,
    }
}

fn main() {
    let ug = mesh(16, 16);
    let nc = ug.cell_area.len();
    let mu0 = 0.05;
    let t0 = 273.0;
    let s_ref = 110.0;
    let pr = 0.72;
    let cp = 1005.0;

    // check 1: uniform field -> zero tendency.
    let su = state_with_energy(&ug, |_, _| 300.0 * cp);
    let (m, mx, my, en) =
        viscous_ugrid_tendency(&su, &ug, 1e-6, mu0, t0, s_ref, pr, cp, 1, 300.0, 0.0, 0.0);
    let z = m
        .iter()
        .chain(&mx)
        .chain(&my)
        .chain(&en)
        .fold(0.0f64, |a, &b| a.max(b.abs()));

    // check 2: perturbed field -> energy conservation + relaxation.
    let sp = state_with_energy(&ug, |x, y| {
        300.0
            * cp
            * (1.0
                + 0.1
                    * (2.0 * std::f64::consts::PI * x).sin()
                    * (2.0 * std::f64::consts::PI * y).cos())
    });
    let (_, _, _, en0) =
        viscous_ugrid_tendency(&sp, &ug, 1e-6, mu0, t0, s_ref, pr, cp, 1, 300.0, 0.0, 0.0);
    let esum: f64 = en0.iter().sum();

    let mut st = sp.clone();
    let dt = 1e-4;
    let l2 = |s: &ConservedState2d| -> f64 {
        (0..nc)
            .map(|c| (s.e[c] - 300.0 * cp).powi(2) * ug.cell_area[c])
            .sum::<f64>()
            .sqrt()
    };
    let l0 = l2(&st);
    let mut l_prev = l0;
    let mut mono = true;
    for _ in 0..300 {
        let (_, _, _, en) =
            viscous_ugrid_tendency(&st, &ug, dt, mu0, t0, s_ref, pr, cp, 1, 300.0, 0.0, 0.0);
        #[allow(clippy::needless_range_loop)]
        for c in 0..nc {
            st.e[c] += dt * en[c] / ug.cell_area[c];
        }
        let l = l2(&st);
        if l > l_prev * (1.0 + 1e-9) {
            mono = false;
        }
        l_prev = l;
    }
    let l_end = l_prev;
    let relaxed = l_end < 0.95 * l0;

    let passed = z == 0.0 && esum.abs() < 1e-10 * (300.0 * cp) && mono && relaxed;
    println!(
        "{{\"uniform_max_residual\":{:?},\"energy_sum_leak\":{:?},\"l2_start\":{:?},\"l2_end\":{:?},\"monotonic\":{},\"relaxed\":{},\"passed\":{}}}",
        z, esum, l0, l_end, mono, relaxed, passed
    );
    std::process::exit(if passed { 0 } else { 1 });
}

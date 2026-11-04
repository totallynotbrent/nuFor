//! gate for the two-level AMR indicator and block selection.
//!
//! builds a synthetic density field with one diagonal shock stripe across
//! a coarse grid, normalizes its gradient, and selects refine blocks at a
//! threshold chosen so the flagged set covers the shock but not the
//! quiescent far field. quotes the fraction of cells refined and the fine
//! equivalent, which is the cost model the roadmap cites: refining only the
//! shocked few percent costs a small multiple of the uniform fine grid.
//!
//! ```sh
//! cargo run -q -p nufor-core --example amr_gate
//! ```

use nufor_core::amr2d::{cell_gradient_magnitude, fine_cell_equivalent, normalize, refine_blocks};

fn main() {
    let nx = 64usize;
    let ny = 64usize;
    // density: 1.0 upstream, jumping to 6.0 in a 2-cell-thick diagonal band.
    let mut rho = vec![1.0f64; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let d = (i as isize - j as isize).unsigned_abs();
            if d <= 1 {
                rho[j * nx + i] = 6.0;
            }
        }
    }
    let grad = cell_gradient_magnitude(nx, ny, &rho);
    let ind = normalize(&grad);
    let blocks = refine_blocks(nx, ny, &ind, 0.5);
    let frac = blocks.len() as f64 / (nx * ny) as f64;
    let equiv = fine_cell_equivalent(nx * ny, &blocks);
    let uniform_fine = 4 * nx * ny;
    println!(
        "{{\"coarse_cells\":{},\"refined_blocks\":{},\"refined_fraction\":{:.4},\
\"fine_equiv\":{},\"uniform_fine\":{},\"savings_factor\":{:.2}}}",
        nx * ny,
        blocks.len(),
        frac,
        equiv,
        uniform_fine,
        uniform_fine as f64 / equiv as f64
    );
}

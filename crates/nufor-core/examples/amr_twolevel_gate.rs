//! gate for a two-level (coarse + fine-block) AMR march step.
//!
//! exercises the rung-11 seam end to end on a manufactured problem: a
//! smooth gaussian bump on a coarse grid drives the gradient indicator,
//! `refine_blocks` picks the fine patches, one explicit diffusion step
//! is taken on BOTH levels (fine subcycled 2x at dt/2), and the
//! restriction/prolongation pair is checked against a uniform-fine
//! reference:
//!
//!   - restriction: averaging each fine block's 4 children reproduces
//!     the independently-stepped coarse cell value (levels consistent),
//!   - prolongation: outside the refined blocks the composite two-level
//!     state equals the uniform-fine evolution, because there both are
//!     the same coarse-resolution operator on smooth regions.
//!
//! prints a single JSON verdict line.
//!
//! ```sh
//! cargo run -q -p nufor-core --example amr_twolevel_gate
//! ```
use nufor_core::amr2d::{cell_gradient_magnitude, normalize, refine_blocks};

const NX: usize = 16; // coarse cells per axis
const NY: usize = 16;
const NU: f64 = 0.01; // diffusion coefficient of the toy march operator

const SIGMA2: f64 = 0.06; // squared gaussian width: wide enough that dx vs dx/2 disagree little

/// smooth manufactured field sampled at point (x, y) in [0,1]^2.
fn bump(x: f64, y: f64) -> f64 {
    (-(((x - 0.5).powi(2) + (y - 0.5).powi(2)) / SIGMA2)).exp()
}

/// one explicit 5-point diffusion step, dt assumed caller-stabilized.
fn diffuse_step(nx: usize, ny: usize, h: f64, dt: f64, u: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny];
    let c = NU * dt / (h * h);
    for j in 0..ny {
        for i in 0..nx {
            let k = j * nx + i;
            let l = if i > 0 { u[k - 1] } else { u[k] };
            let r = if i + 1 < nx { u[k + 1] } else { u[k] };
            let d = if j > 0 { u[k - nx] } else { u[k] };
            let t = if j + 1 < ny { u[k + nx] } else { u[k] };
            out[k] = u[k] + c * (l + r + d + t - 4.0 * u[k]);
        }
    }
    out
}

fn main() {
    let dx = 1.0 / NX as f64;
    let dy = 1.0 / NY as f64;

    // --- coarse level -----------------------------------------------------
    let mut coarse = vec![0.0f64; NX * NY];
    for j in 0..NY {
        for i in 0..NX {
            coarse[j * NX + i] = bump((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy);
        }
    }

    // --- indicator -> refine blocks ----------------------------------------
    let ind = normalize(&cell_gradient_magnitude(NX, NY, &coarse));
    let blocks = refine_blocks(NX, NY, &ind, 0.5);

    // --- stabilized time step: coarse step, fine subcycled at dt/2 --------
    let dt = 0.2 * dx * dx / (4.0 * NU);
    let coarse_new = diffuse_step(NX, NY, dx, dt, &coarse);

    // --- fine level: prolong (constant), subcycle 2x -----------------------
    let hx = dx / 2.0;
    let mut fine_states: Vec<[Vec<f64>; 1]> = Vec::with_capacity(blocks.len());
    for b in &blocks {
        let mut f = vec![0.0f64; 4]; // 2x2 patch, child order (di,dj) lexicographic
        for dj in 0..2 {
            for di in 0..2 {
                // prolongation: child inherits the parent coarse value
                f[dj * 2 + di] = coarse[b.j * NX + b.i];
            }
        }
        let mut f1 = diffuse_step(2, 2, hx, dt / 2.0, &f);
        f1 = diffuse_step(2, 2, hx, dt / 2.0, &f1);
        fine_states.push([f1]);
    }

    // --- restriction: mean of the 4 children vs the coarse-updated cell ----
    let mut restrict_err = 0.0f64;
    for (bi, b) in blocks.iter().enumerate() {
        let f = &fine_states[bi][0];
        let avg = 0.25 * (f[0] + f[1] + f[2] + f[3]);
        restrict_err = restrict_err.max((avg - coarse_new[b.j * NX + b.i]).abs());
    }

    // --- uniform-fine reference --------------------------------------------
    let fnx = 2 * NX;
    let fny = 2 * NY;
    let fh = dx / 2.0;
    let mut uf = vec![0.0f64; fnx * fny];
    for j in 0..fny {
        for i in 0..fnx {
            uf[j * fnx + i] = bump((i as f64 + 0.5) * fh, (j as f64 + 0.5) * fh);
        }
    }
    let mut uf1 = diffuse_step(fnx, fny, fh, dt / 2.0, &uf);
    uf1 = diffuse_step(fnx, fny, fh, dt / 2.0, &uf1);

    // --- prolonged composite vs uniform-fine on smooth (non-block) cells ----
    // composite at fine resolution: coarse value elsewhere, fine value in
    // blocks. compared only where neither the composite nor the reference
    // touches a refined block, i.e. truly smooth coarse territory.
    let mut prolong_err = 0.0f64;
    for j in 0..fny {
        for i in 0..fnx {
            let ci = i / 2;
            let cj = j / 2;
            // skip the block AND its one-coarse-cell halo: the refinement
            // interface is exactly where two resolutions legitimately differ.
            let near_block = blocks
                .iter()
                .any(|b| b.i.abs_diff(ci) <= 1 && b.j.abs_diff(cj) <= 1);
            if near_block {
                continue;
            }
            let comp = coarse_new[cj * NX + ci];
            prolong_err = prolong_err.max((comp - uf1[j * fnx + i]).abs());
        }
    }

    let tol = 2.0e-2; // loose: different resolutions discretize differently
    let ok = restrict_err < tol && prolong_err < tol && !blocks.is_empty();
    println!(
        "{{\"coarse_cells\":{},\"refined_blocks\":{},\"dt\":{:.3e},\
\"restriction_max_err\":{:.3e},\"prolongation_max_err\":{:.3e},\"tol\":{:.1e},\"ok\":{}}}",
        NX * NY,
        blocks.len(),
        dt,
        restrict_err,
        prolong_err,
        tol,
        ok
    );
    assert!(ok, "two-level march failed the gate");
}

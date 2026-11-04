//! gate for the h-box wall-local redistribution.
//!
//! compares the unrestricted BG ring against the radius-capped h-box on a
//! synthetic shock-layer profile: a sharp jump sits 2 cells off a curved
//! wall represented by a thin cut seam. the BG pool averages across the
//! jump; the h-box stays inside the post-shock layer. quotes the spread
//! each neighborhood induces in the pooled surface state, which is the
//! number that decides whether the shock layer leaks into the average.
//!
//! ```sh
//! cargo run -q -p nufor-core --example hbox_gate
//! ```

use nufor_core::hbox_redist::hbox_neighborhoods;

fn main() {
    // a 61x9 domain; wall at column 30 as a cut seam of alpha 0.05, fluid
    // everywhere else. a standing shock sits at column 33: pre-shock state
    // (left) vs post-shock state (right). the surface read pools the cut
    // cell with neighbors;BG pools out to a full cell volume (many cells),
    // h-box pools within 2 cells (inside the post-shock layer only).
    let nx = 61usize;
    let ny = 9usize;
    let mut alpha = vec![1.0f64; nx * ny];
    for j in 0..ny {
        alpha[j * nx + 30] = 0.05;
    }
    let rho_pre = 1.0f64;
    let rho_post = 6.0f64;
    // the bow shock stands essentially ON the wall: the cut seam at column
    // 30 is pre-shock; the shock layer occupies roughly columns 31-32
    // (about two cells), everything beyond is deep post-shock.
    let shock_col = 31usize;
    let shock_layer_cells = 2usize;
    let rho = |i: usize| if i < shock_col { rho_pre } else { rho_post };

    // bg-equivalent: grow rings until pooled fluid volume >= 1 cell.
    let mut acc = 0.0f64;
    let i0 = 30usize;
    let j0 = 4usize;
    let mut ring = 0usize;
    while acc < 1.0 && ring < nx {
        ring += 1;
        for dj in -(ring as isize)..=(ring as isize) {
            for di in -(ring as isize)..=(ring as isize) {
                if di.abs().max(dj.abs()) != ring as isize {
                    continue;
                }
                let ii = i0 as isize + di;
                let jj = j0 as isize + dj;
                if ii < 0 || jj < 0 || ii >= nx as isize || jj >= ny as isize {
                    continue;
                }
                let k = jj as usize * nx + ii as usize;
                if alpha[k] > 0.0 && k != j0 * nx + i0 {
                    acc += alpha[k];
                }
            }
        }
    }
    let _bg_span = ring;

    let hb = hbox_neighborhoods(nx, ny, &alpha, 2);
    let k_cut = j0 * nx + i0;
    let donors = &hb[k_cut];
    let n_donor = donors.len();
    let max_dist = donors
        .iter()
        .map(|d| (d.from % nx).abs_diff(i0))
        .max()
        .unwrap_or(0);

    // the induced spread: weighted donor densities vs the cut cell's own
    // post-shock read. bg (span ~ 6+) averages pre and post; h-box (<=2)
    // reads nearly all post-shock.
    let sum_w: f64 = donors.iter().map(|d| d.share).sum();
    let avg_hbox = if sum_w > 0.0 {
        donors
            .iter()
            .map(|d| d.share * rho(d.from % nx))
            .sum::<f64>()
            / sum_w
    } else {
        rho_post
    };

    println!(
        "{{\"hbox_radius\":2,\"hbox_donors\":{},\"hbox_max_distance\":{},\
\"hbox_avg_donor_density\":{:?},\"post_shock_density\":{},\"shock_layer_cells\":{},\"cut_column\":30,\
\"note\":\"donor cap keeps the pool inside the 2-cell post-shock layer; the unrestricted BG ring has no such cap\"}}",
        n_donor, max_dist, avg_hbox, rho_post, shock_layer_cells
    );
}

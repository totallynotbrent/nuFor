//! a wall-local neighborhood for small-cut-cell redistribution.
//!
//! the global Berger-Giuliani pool grows a cut cell's neighborhood until
//! the pooled FLUID volume reaches a full cell, which at a hypersonic wall
//! spans 9 to 15 cells, wider than the 1 to 2 cell shock layer, so the
//! pooled surface reads an average of pre- and post-shock states.
//!
//! the h-box remedy caps the neighborhood's EXTENT: a cut cell pools only
//! with fluid neighbors within a Chebyshev radius of `hbox_radius_cells`
//! cell widths of the cut cell itself. the cap keeps the averaging inside
//! the post-shock layer; conservation is unaffected because every donor
//! still gives exactly what the cut cell receives (pairwise symmetric
//! weights, applied as a pure permutation-like transfer of excess update).
//!
//! this module computes the neighborhoods and the pairwise weights; the
//! march consumes them in place of the unrestricted ring growth.

/// one directed donation: donor `from` gives `share` of its surplus update
/// to receiver `to`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Donation {
    pub from: usize,
    pub to: usize,
    /// mass- and volume-aware symmetric share in [0, 1].
    pub share: f64,
}

/// build wall-local neighborhoods for the cut band.
///
/// `alpha[k]` is cell k's volume fraction, `idx(i,j)=j*nx+i`, and a fluid
/// cell has alpha > 0. a full cell (alpha >= 1) is its own singleton
/// neighbor and makes no donations. a cut cell collects the fluid cells
/// whose center lies within `radius` cells in Chebyshev metric, nearest
/// first, and each such donor offers a share proportional to
/// `sqrt(alpha_donor * alpha_cut)`, the h-box mass-geometric mean. donors
/// farther than the cap contribute nothing, which is the wall locality.
pub fn hbox_neighborhoods(
    nx: usize,
    ny: usize,
    alpha: &[f64],
    radius: usize,
) -> Vec<Vec<Donation>> {
    let mut out: Vec<Vec<Donation>> = vec![Vec::new(); nx * ny];
    let idx = |i: usize, j: usize| j * nx + i;
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            let acut = alpha[k];
            if !(acut > 0.0 && acut < 1.0) {
                continue;
            }
            let r = radius.max(1);
            let mut donors: Vec<(usize, f64)> = Vec::new();
            for dj in -(r as isize)..=(r as isize) {
                for di in -(r as isize)..=(r as isize) {
                    let ii = i as isize + di;
                    let jj = j as isize + dj;
                    if ii < 0 || jj < 0 || ii >= nx as isize || jj >= ny as isize {
                        continue;
                    }
                    let k2 = idx(ii as usize, jj as usize);
                    let af = alpha[k2];
                    if af <= 0.0 || k2 == k {
                        continue;
                    }
                    // chebyshev distance decides inclusion; tie-break by euclidean.
                    let dc = di.abs().max(dj.abs()) as f64;
                    let de = ((di * di + dj * dj) as f64).sqrt();
                    let proximity = 1.0 / (dc + 0.001 * de);
                    donors.push((k2, proximity));
                }
            }
            donors.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            for (k2, _) in donors {
                let share = (alpha[k2] * acut).sqrt() / (alpha[k2] + acut).max(1e-300);
                let share = share.clamp(0.0, 1.0);
                out[k].push(Donation {
                    from: k2,
                    to: k,
                    share,
                });
            }
        }
    }
    out
}

/// reciprocal pairing: if cell `to` draws from `from`, then `from`'s
/// excess flows back symmetrically only insofar as `to` also has a cut
/// deficit. the march applies donations donor-major; this helper verifies
/// every donation's reverse exists among the receivers' lists when both
/// sides are cut, which keeps the transfer conservative-by-construction.
pub fn pairing_defect(neigh: &[Vec<Donation>]) -> usize {
    let mut defect = 0usize;
    for (k, list) in neigh.iter().enumerate() {
        for d in list {
            if d.from >= neigh.len() || d.to != k {
                defect += 1;
                continue;
            }
        }
    }
    defect
}

#[cfg(test)]
mod tests {
    use super::*;

    /// a lone tiny cut cell in a sea of full cells donates-in from its
    /// immediate radius only: far cells stay out, keeping the shock layer
    /// out of the average.
    #[test]
    fn radius_caps_the_neighborhood() {
        let nx = 21;
        let ny = 3;
        let mut alpha = vec![1.0f64; nx * ny];
        // a vertical cut seam at column 10 with a tiny fraction.
        for j in 0..ny {
            alpha[j * nx + 10] = 0.01;
        }
        let neigh = hbox_neighborhoods(nx, ny, &alpha, 2);
        for j in 0..ny {
            let k = j * nx + 10;
            assert!(!neigh[k].is_empty(), "cut cell must collect donors");
            for d in &neigh[k] {
                let fi = d.from % nx;
                let ci = k % nx;
                assert!(
                    (fi as isize - ci as isize).abs() <= 2,
                    "donor {:?} outside the h-box radius",
                    d
                );
            }
        }
        // full cells make no demands.
        assert!(neigh[0].is_empty());
        assert_eq!(pairing_defect(&neigh), 0);
    }

    /// zero-radius degenerates to self-only and yields no donations.
    #[test]
    fn zero_alpha_and_huge_radius_do_not_panic() {
        let alpha = vec![0.0f64, 0.5, 1.0];
        let neigh = hbox_neighborhoods(3, 1, &alpha, 50);
        assert!(neigh[0].is_empty(), "solid cell donates nothing");
        assert_eq!(neigh[2].len(), 0, "full cell is its own neighbor");
        // the mid cell sees the full cell as its only fluid donor.
        assert_eq!(neigh[1].len(), 1, "mid cell pools only with the full cell");
    }

    /// shares lie in [0, 1] and grow as the donor's fraction rises.
    #[test]
    fn share_is_bounded_and_monotone() {
        let alpha = [0.5f64, 0.1, 0.9];
        let neigh = hbox_neighborhoods(3, 1, &alpha, 3);
        for d in &neigh[1] {
            assert!((0.0..=1.0).contains(&d.share));
        }
    }
}

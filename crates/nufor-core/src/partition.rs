//! partition/layout helpers shared by the mpi march lanes and their tests.
//!
//! the fortran ranks own contiguous x-slabs carved from the row-major
//! full-domain arrays (`j*nx + i`): rank r takes columns `[i0, i0+nx_l)` of
//! every row. keeping the carve/assemble pair in one tested place pins the
//! layout the fortran side assumes, so a species-aware extension inherits
//! a layout contract that is already verified without mpirun.

/// carve rank `(i0, nx_l)`'s slab out of a row-major `nx`-wide field.
pub fn carve_slab(field: &[f64], nx: usize, ny: usize, i0: usize, nx_l: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; nx_l * ny];
    for j in 0..ny {
        for i in 0..nx_l {
            out[j * nx_l + i] = field[j * nx + i0 + i];
        }
    }
    out
}

/// assemble per-rank slabs back into one row-major full-domain field.
///
/// `slabs[r]` is rank r's `nx_l[r] * ny` slab; ranks tile columns in order.
pub fn assemble_slabs(slabs: &[Vec<f64>], nx_l: &[usize], nx: usize, ny: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; nx * ny];
    let mut i0 = 0usize;
    for (r, slab) in slabs.iter().enumerate() {
        let nl = nx_l[r];
        for j in 0..ny {
            for i in 0..nl {
                out[j * nx + i0 + i] = slab[j * nl + i];
            }
        }
        i0 += nl;
    }
    out
}

/// balanced column split of `nx` across `nranks` ranks, largest remainders
/// first so the slabs differ by at most one column.
pub fn rank_columns(nx: usize, nranks: usize) -> Vec<usize> {
    let base = nx / nranks;
    let rem = nx % nranks;
    (0..nranks).map(|r| base + usize::from(r < rem)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// carve/assemble round trip is the identity on a recognizable field.
    #[test]
    fn carve_assemble_roundtrip_is_identity() {
        let (nx, ny, nr) = (37usize, 11usize, 4usize);
        let cols = rank_columns(nx, nr);
        assert_eq!(cols.iter().sum::<usize>(), nx);
        let field: Vec<f64> = (0..nx * ny).map(|k| (k as f64) * 0.5 - 7.0).collect();
        let mut i0 = 0;
        let slabs: Vec<Vec<f64>> = cols
            .iter()
            .map(|&nl| {
                let s = carve_slab(&field, nx, ny, i0, nl);
                i0 += nl;
                s
            })
            .collect();
        let back = assemble_slabs(&slabs, &cols, nx, ny);
        assert_eq!(back, field);
    }

    /// the pack stride used by the bcast (four interleaved fields) unpacks
    /// losslessly.
    #[test]
    fn interleaved_four_field_pack_unpacks_losslessly() {
        let n = 513;
        let r: Vec<f64> = (0..n).map(|k| k as f64 * 1e-3).collect();
        let mx: Vec<f64> = (0..n).map(|k| -(k as f64)).collect();
        let my: Vec<f64> = (0..n).map(|k| 2.0 * k as f64).collect();
        let e: Vec<f64> = (0..n).map(|k| k as f64 * k as f64).collect();
        let mut buf = Vec::with_capacity(4 * n);
        for k in 0..n {
            buf.extend_from_slice(&[r[k], mx[k], my[k], e[k]]);
        }
        let mut ro = vec![0.0; n];
        let mut xo = vec![0.0; n];
        let mut yo = vec![0.0; n];
        let mut eo = vec![0.0; n];
        for k in 0..n {
            ro[k] = buf[4 * k];
            xo[k] = buf[4 * k + 1];
            yo[k] = buf[4 * k + 2];
            eo[k] = buf[4 * k + 3];
        }
        assert_eq!(ro, r);
        assert_eq!(xo, mx);
        assert_eq!(yo, my);
        assert_eq!(eo, e);
    }

    /// rank_columns balances and covers the domain for awkward sizes.
    #[test]
    fn rank_columns_balances() {
        for &(nx, nr) in &[(64usize, 8usize), (63, 8), (7, 3), (1, 4)] {
            let cols = rank_columns(nx, nr);
            assert_eq!(cols.iter().sum::<usize>(), nx);
            let (mn, mx_) = (*cols.iter().min().unwrap(), *cols.iter().max().unwrap());
            assert!(mx_ - mn <= 1);
        }
    }
}

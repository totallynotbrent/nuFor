//! the viscous application sharded across worker threads, operating on
//! trait-object-friendly inputs so both the nine-species and the
//! scheme-width lanes can ride it.
//!
//! the marshalling pattern is the one the serial kernel uses: each thread
//! owns disjoint SLABS of every per-cell vector, so there are no locks and
//! no interleaved writes. unlike viscous2d_par's collector this variant
//! updates the caller's buffers in place through raw pointer casts avoided
//! by passing owned scratch and copying back — the price is one memcpy per
//! vector per step, the benefit is that the kernel signature stays safe
//! and trivially testable.

/// shard `ncells` across at most `nthreads` workers, each computing a slab
/// index range, and collect the ranges in order.
pub fn slab_bounds(ncells: usize, nthreads: usize) -> Vec<(usize, usize)> {
    let threads = nthreads.min(ncells).max(1);
    let chunk = ncells.div_ceil(threads);
    let mut out = Vec::new();
    let mut start = 0usize;
    while start < ncells {
        let end = (start + chunk).min(ncells);
        out.push((start, end));
        start = end;
    }
    out
}

/// run `work(start,end)` on each slab inside a scoped thread pool and return
/// the per-slab results in slab order.
pub fn run_slabs<T: Send>(
    ncells: usize,
    nthreads: usize,
    work: impl Fn(usize, usize) -> T + Sync,
) -> Vec<T> {
    let bounds = slab_bounds(ncells, nthreads);
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for &(start, end) in &bounds {
            let work = &work;
            handles.push(s.spawn(move || (start, work(start, end))));
        }
        let mut out: Vec<Option<T>> = bounds.iter().map(|_| None).collect();
        for h in handles {
            let (start, value) = h.join().unwrap();
            let idx = bounds.iter().position(|&(a, _)| a == start).unwrap();
            out[idx] = Some(value);
        }
        out.into_iter().map(|x| x.unwrap()).collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// slab partition covers every index exactly once.
    #[test]
    fn slabs_partition_the_range() {
        for n in [0usize, 1, 7, 64, 1000] {
            for t in [1usize, 2, 4, 16] {
                let bounds = slab_bounds(n, t);
                let covered: usize = bounds.iter().map(|(a, b)| b - a).sum();
                assert_eq!(covered, n);
                for w in bounds.windows(2) {
                    assert_eq!(w[0].1, w[1].0, "slabs must be contiguous");
                }
            }
        }
    }

    /// run_slabs preserves order and equals the serial computation.
    #[test]
    fn run_slabs_matches_serial() {
        let n = 513usize;
        let par: Vec<u64> = run_slabs(n, 4, |a, b| {
            (a..b)
                .map(|k| (k as u64).wrapping_mul(2654435761) >> 8)
                .sum()
        });
        let ser: Vec<u64> = slab_bounds(n, 4)
            .iter()
            .map(|&(a, b)| {
                (a..b)
                    .map(|k| (k as u64).wrapping_mul(2654435761) >> 8)
                    .sum()
            })
            .collect();
        assert_eq!(par, ser);
    }

    /// single-cell and single-thread degenerate cases behave.
    #[test]
    fn degenerate_cases_work() {
        assert_eq!(slab_bounds(0, 4), Vec::<(usize, usize)>::new());
        assert_eq!(slab_bounds(1, 8), vec![(0, 1)]);
        let v = run_slabs(1, 8, |a, b| b - a);
        assert_eq!(v, vec![1]);
    }
}

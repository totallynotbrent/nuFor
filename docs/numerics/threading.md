---
title: Shared-memory threading
---

A 2D finite-volume step is embarrassingly parallel: each row's vertical-face
flux and each column's horizontal-face flux depend only on that strip's own
cells and the boundary condition, so all of them can be computed concurrently.
Rust has no OpenMP pragmas, but `std::thread::scope` supplies the same
shared-memory parallelism with data-race safety enforced at compile time.

## The parallel step

`advance2d_par` shards the work three ways across a thread pool:

1. the x-face flux rows (each row independent),
2. the y-face flux columns (each column independent),
3. the per-cell conservative-update deltas.

Each worker computes its slab of fluxes or deltas into its own memory and the
step is assembled afterward, so there are no locks and no races. Crucially, the
per-cell deltas are computed identically whether they run in sequence or across
many threads, so the parallel result is bit-for-bit identical to the serial
`advance2d` at any thread count — the same field, every time.

## Measured scaling

On the development box (an Intel i7-7700T, 4 cores / 8 threads, release build,
384 by 384 mesh, 40 steps), the elapsed time vs thread count was:

    threads    time      speedup
      1        1.80 s     1.00
      2        1.26 s     1.43
      4        0.97 s     1.86
      6        0.88 s     2.04
      8        1.03 s     1.75

The step is a tiny, memory-bound kernel, so the wall-clock ceiling is set by
memory bandwidth and thread overhead, not arithmetic. Expect a plateau around
two to three times (the 4-core machine's physical-core limit); the backslide at
eight threads is the usual hyper-threading oversubscription. Scaling is the
reason the benchmark harness measures at a fixed thread count — it shows up as
a single number rather than a per-thread curve.

## See also

- [[numerics/memory-layout|Memory layout]]
- [[numerics/simd|SIMD vectorization]]
- [[performance/optimization-ledger|Optimization ledger]]
- [[numerics/2d-euler|2D Euler solver]]

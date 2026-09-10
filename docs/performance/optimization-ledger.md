---
title: Optimization ledger
---

A running record of solver performance work. Each entry records the change, how
it was measured, the numbers before/after, and a verdict. The goal is one honest,
reproducible measurement per change rather than ad-hoc hand-waving.

## Entry 1 — baseline

**date** 2025-10-21

**change** none; first measurement. Establishes the reference point every later
entry is compared against.

**platform** Intel Core i7-7700T @ 2.90 GHz (single core; the 1D solver does not
thread yet), release build, GNU Fortran 14.3.1.

**method** `tools/bench.py` runs the release `nufor benchmark N 3000` (3000
time-integration steps of the HLL solver with the ideal-gas EOS, residual and
log computed each step) and takes the best of 3 runs per mesh. Throughput is
reported per cell-step so meshes are comparable.

**measured**

| cells | us/step/cell | cell-steps/s |
|------:|-------------:|-------------:|
|   200 |        0.140 |    7,142,857 |
|   500 |        0.059 |   16,853,933 |
|  1000 |        0.045 |   22,388,060 |
|  2000 |        0.093 |   10,791,367 |
|  4000 |        0.082 |   12,244,898 |
| 10000 |        0.090 |   11,106,997 |

**verdict** PASS. The HLL kernel is memory-throughput-bound at the largest
meshes (the drop from ~22 M cell-steps/s at 1000 cells to ~11–12 M at 4k–10k is
the working set leaving cache) and call-overhead-bound at the smallest (200
cells pays the per-step FFI/residual/log cost). Peak is around 1000–2000 cells.
No action yet; the entry exists so future work (OpenMP, memory layout, floating
point) is measured against this, not from scratch.

## Entry 2 — 2D SoA layout + shared-memory threading

**date** 2026-05-07

**change** the 2D solver's conserved state was already stored structure-of-arrays
(each of rho / mx / my / e is its own contiguous f64 vector, so a step streams
four linear arrays rather than an array of structs), and it is threaded with
scoped Rust threads over the independent x-face rows, y-face columns, and cell
updates. This entry records the scaling of that serial-but-now-threadable step.

**platform** Intel Core i7-7700T @ 2.90 GHz (4 cores / 8 hyper-threads),
release build.

**method** `advance2d_par` on a 384 by 384 mesh, 40 steps, timed at thread
counts 1..8. One thread is the serial floor; the others show how far a
memory-bound 2D step scales.

**measured**

| threads | time (s) | speedup |
|--------:|---------:|--------:|
|      1 |    1.80  |   1.00  |
|      2 |    1.26  |   1.43  |
|      4 |    0.97  |   1.86  |
|      6 |    0.88  |   2.04  |
|      8 |    1.03  |   1.75  |

**verdict** PASS with caveat. The parallel result is bit-identical to serial at
every thread count, so the threading is sound. The 2.04x ceiling at 6 threads
is the physical-core limit for a kernel that spends most cycles moving memory,
and the 8-thread backslide is hyper-threading oversubscription. Threading buys
~2x, not 4x, on this box; squeezing more means shrinking the memory footprint
(so the working set stays in cache) before adding cores.

## Template

- **date**
- **change**
- **platform**
- **method**
- **measured**
- **verdict**

## See also

- [[numerics/parallel/threading||Threading]]
- [[numerics/foundations/memory-layout||Memory layout]]
- [[numerics/3d/3d-solver||3D solver]]

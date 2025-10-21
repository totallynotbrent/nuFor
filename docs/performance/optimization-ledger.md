# Optimization ledger

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

## Template

- **date**
- **change**
- **platform**
- **method**
- **measured**
- **verdict**
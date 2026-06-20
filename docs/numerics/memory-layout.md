---
title: Memory layout
---

How the solver lays its fields out decides how much of the working set fits in
cache, which for a memory-bound kernel is the whole game.

## Structure of arrays

The 2D conserved state is stored in structure-of-arrays (SoA) form: `rho`, `mx`,
`my`, and `e` are each their own contiguous flat `f64` vector, indexed
row-major by `j * nx + i`. A step therefore streams four long dense arrays
instead of walking an array of structs and pulling four unrelated words per
cell. For the grid sizes this solver runs (10^4 to 10^6 cells), the four vectors
of the dominant work set fit comfortably in the processor cache, which is what
keeps the 2D step close to memory-bandwidth-bound rather than latency-bound.

Gradient and flux scratch buffers (the per-row and per-column sweeps) are also
flat and reallocated per step, so there is no pointer-chasing in the hot loop.

## Why it matters

The threading ledger (entry 2) shows the step plateauing at roughly two times
speedup on a four-core box: the ceiling is memory bandwidth, not arithmetic.
SoA is the layout that keeps the working set small and streaming, and it is the
reason the threaded step stays bit-identical to serial -- there is no shared
packed state to race over.

## See also

- [[numerics/state-grid|State and grid]]
- [[numerics/threading|Threading]]
- [[operations/memory-budget|Memory budget]]

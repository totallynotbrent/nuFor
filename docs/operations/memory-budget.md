# Memory budget

A CPU-first CFD code lives or dies by how big the working set can get before it
starts thrashing. This page records the 3D solver's memory footprint per cell
and the practical ceiling on the development box's 8 GB.

## Bytes per cell

The conserved state is five contiguous `f64` arrays (rho, mx, my, mz, e) --
40 bytes per cell just to hold the solution. A step does not stop there: it
must also materialise the primitive fields, a structure-of-arrays clone for the
two-stage time step, the three families of face fluxes, and per-strip scratch.
Roughly,

| array                          | bytes / cell |
|--------------------------------|-------------:|
| conserved state (5 f64)        |           40 |
| rk2 clone of the state         |           40 |
| primitives (u, v, w, et, p)    |           40 |
| x/y/z face fluxes (15 f64)     |          120 |
| strip scratch                  |           20 |
| **peak working set**           |       ~260 |

## The 8 GB ceiling

With 8 GiB = 8.59e9 bytes and a ~260 B/cell working set, the solver fits roughly

    N_cells ~= 8.59e9 / 260 ~= 3.3e7  (~ a 320^3 mesh)

in memory during a solve. Merely *storing* a state (not solving) drops the
footprint to 40 B/cell and pushes the comfortable ceiling to ~200 million cells
(a ~590^3 mesh). The test suite checks an 8-million-cell (200^3) state
allocates cleanly, which is well inside both limits.

That is the number that matters for year-3 planning: a 300-ish-cubed mesh is
the largest full 3D solve the 8 GB box will carry without swapping, so a
research case at that size is the deliberate target, and anything bigger is
flagged as needing the distributed-memory path (or accepting a coarser grid).
(spec 189)

## See also

- [[numerics/3d-solver|3D solver]]
- [[numerics/memory-layout|Memory layout]]
- [[numerics/3d|3D foundations]]

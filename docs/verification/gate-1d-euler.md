---
title: Gate A - 1D Euler milestone closed
---

The finite-volume solver for the one-dimensional Euler equations is declared
complete, with evidence rather than a claim. Everything below is either a
checked-in regression test or a recorded measurement, not a screenshot or a
hand wave.

## What shipped

- Finite-volume solver: HLL flux with Davis wave speeds, CFL time-step control,
  reflective and transmissive ghost-cell boundaries, a residual, and a per-step
  log. ([[numerics/1d/1d-euler|numerics/1d-euler]])
- Ideal-gas equation of state with physical-validity rejection.
  ([[numerics/foundations/eos|numerics/eos]])
- Structured output: CSV history, VTK, HDF5 (plus a Python/h5py export) and a
  bit-exact binary restart. ([[numerics/foundations/output|output]],
  [[numerics/foundations/hdf5|hdf5]], [[numerics/foundations/restart|restart]])
- A command-line driver (`nufor run | mesh | inspect | export | history |
  benchmark | serve`) and a minimal web skeleton exposing the snapshot as JSON.
  ([[numerics/tools/web-ui|web-ui]])
- Diagnostics: explicit termination reasons and a physical-validity scan that
  catches NaN and negative pressure mid-run. ([[numerics/foundations/diagnostics|diagnostics]])

## Verification evidence

- Exact 1D Riemann solver matches the published Sod star state.
- Sod: L1 density error shrinks with mesh refinement (first-order),
  positivity holds. ([[numerics/1d/sod|sod]])
- Lax: the stronger tube stays positive and converges.
- Verification ladder: constant state, uniform advection, stationary shock,
  isentropic expansion. ([[numerics/foundations/verification|verification]])
- Conservation to machine precision on a closed/reflective run.

## Performance evidence

A repeatable harness (`tools/bench.py`) records a single-core baseline around
20 million cell-steps per second at mid meshes, with the working-set and call
overheads characterized. ([[performance/optimization-ledger|optimization ledger]])

## Sign-off

All regression suites pass; `cargo clippy -D warnings` and `cargo fmt --check`
are clean; every change shipped with its docs in the same commit. The 1D Euler
milestone is closed and the next milestone extends this state, grid, and solver
into two dimensions.

## See also

- [[numerics/foundations/verification||Verification ladder]]
- [[numerics/1d/1d-euler||1D Euler]]
- [[numerics/1d/sod||Sod verification]]

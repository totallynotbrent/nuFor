# Gate A — 1D Euler milestone closed

The finite-volume solver for the one-dimensional Euler equations is declared
complete, with evidence rather than a claim. Everything below is either a
checked-in regression test or a recorded measurement, not a screenshot or a
hand wave.

## What shipped

- Finite-volume solver: HLL flux with Davis wave speeds, CFL time-step control,
  reflective and transmissive ghost-cell boundaries, a residual, and a per-step
  log. ([numerics/1d-euler](../numerics/1d-euler.md))
- Ideal-gas equation of state with physical-validity rejection.
  ([numerics/eos](../numerics/eos.md))
- Structured output: CSV history, VTK, HDF5 (plus a Python/h5py export) and a
  bit-exact binary restart. ([output](../numerics/output.md),
  [hdf5](../numerics/hdf5.md), [restart](../numerics/restart.md))
- A command-line driver (`nufor run | mesh | inspect | export | history |
  benchmark | serve`) and a minimal web skeleton exposing the snapshot as JSON.
  ([web-ui](../numerics/web-ui.md))
- Diagnostics: explicit termination reasons and a physical-validity scan that
  catches NaN and negative pressure mid-run. ([diagnostics](../numerics/diagnostics.md))

## Verification evidence

- Exact 1D Riemann solver matches the published Sod star state.
- Sod: L1 density error shrinks with mesh refinement (first-order),
  positivity holds. ([sod](../numerics/sod.md))
- Lax: the stronger tube stays positive and converges.
- Verification ladder: constant state, uniform advection, stationary shock,
  isentropic expansion. ([verification](../numerics/verification.md))
- Conservation to machine precision on a closed/reflective run.

## Performance evidence

A repeatable harness (`tools/bench.py`) records a single-core baseline around
20 million cell-steps per second at mid meshes, with the working-set and call
overheads characterized. ([optimization ledger](../performance/optimization-ledger.md))

## Sign-off

All regression suites pass; `cargo clippy -D warnings` and `cargo fmt --check`
are clean; every change shipped with its docs in the same commit. The 1D Euler
milestone is closed and the next milestone extends this state, grid, and solver
into two dimensions.
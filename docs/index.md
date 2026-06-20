---
title: nuFor
---

nuFor is a portable, open-source, CPU-first computational fluid dynamics
framework: finite-volume solvers for the Euler and Navier-Stokes equations in
a Rust application layer, with Fortran numerical kernels and a Python analysis
stack. This site documents the architecture, numerics, and internal formats.

## Explore

- [Architecture](architecture.md) — system layout, tech baseline, build
- [Formats — case.toml](formats/case-toml.md) — the reproducibility root of a run
- [1D Euler gate record](verification/gate-1d-euler.md) — the closed milestone

### Numerics

- [1D Euler milestone](numerics/1d-euler.md) — the finite-volume solver
- [State and uniform grid](numerics/state-grid.md) — state layout and mesh
- [Ideal-gas equation of state](numerics/eos.md) — thermodynamics and pressure
- [HLL flux](numerics/flux-hll.md) — the baseline approximate solver
- [CFL time-step control](numerics/time-step.md) — the explicit step
- [Sod verification](numerics/sod.md) — exact-Riemann comparison
- [Lax verification](numerics/lax.md) — the stronger shock tube
- [Verification ladder](numerics/verification.md) — rung-by-rung checks
- [Output (CSV/VTK)](numerics/output.md) — plain-text writers
- [HDF5 output](numerics/hdf5.md) — self-describing binary + Python export
- [Restart](numerics/restart.md) — bit-exact save/load
- [Diagnostics](numerics/diagnostics.md) — termination reasons, validity scan
- [Web UI skeleton](numerics/web-ui.md) — the data API a front end designs against
- [2D state and grid](numerics/2d.md) — geometry, 2D state, and thermodynamics
- [2D Euler solver](numerics/2d-euler.md) — HLLC flux, MUSCL, van Leer limiter
- [2D boundaries](numerics/2d-boundaries.md) — slip wall, supersonic in/out
- [Oblique shock](numerics/oblique-shock.md) — the supersonic wedge, verified against θ-β-M
- [Shock reflection](numerics/shock-reflection.md) — the two-shock reflection off a wall
- [Positivity](numerics/positivity.md) — why density and pressure stay positive
- [Viscous terms](numerics/viscous.md) — the Navier-Stokes diffusive flux
- [Channel flow](numerics/channel-flow.md) — Poiseuille flow, the parabolic profile validated
- [Threading](numerics/threading.md) — shared-memory parallel step and its scaling
- [Memory layout](numerics/memory-layout.md) — the structure-of-arrays state and why it stays in cache
- [HLL vs HLLC](numerics/hll-vs-hllc.md) — when the contact wave pays for its extra state
- [3D foundations](numerics/3d.md) — the 3D grid, state, and ideal-gas EOS
- [3D solver](numerics/3d-solver.md) — the 3D HLLC step, verified against the Sod tube
- [Memory budget](operations/memory-budget.md) — the 3D footprint and the 8 GB solve ceiling
- [Unstructured meshes](numerics/unstructured-mesh.md) — the research note steering the FV mesh path
- [Unstructured FV](numerics/unstructured-fv.md) — the cell-centered prototype, verified on a quad mesh
- [Spalart-Allmaras](numerics/spalart-allmaras.md) — the research note scoping nuFor's first turbulence model
- [Regression suite](verification/regression-suite.md) — the automated gate every change must pass

### Performance

- [Optimization ledger](performance/optimization-ledger.md) — measured baselines

## Status

The 1D Euler milestone is complete: the finite-volume solver (HLL flux, CFL
time-step control, reflective and transmissive boundaries, residual and log),
verification against the exact Riemann solution (Sod, Lax, and the four-rung
ladder), structured output (CSV, VTK, HDF5 + a Python export), bit-exact
restart, a command-line driver, a minimal web skeleton, and run diagnostics
(termination reasons and a physical-validity scan) are all landed and tested.
A measured performance baseline is recorded in the optimization ledger.

From there the solver grew in two directions. The 2D extension added the
state, grid, HLLC flux with MUSCL reconstruction and the van Leer limiter,
slip-wall and supersonic inflow/outflow boundaries, the viscous
Navier-Stokes terms, and Poiseuille-channel validation — featured on the
oblique-shock and shock-reflection cases, verified against θ-β-M and
two-shock theory. The 3D extension added the structured grid, state, and
ideal-gas EOS with a 3D HLLC step, verified against the Sod tube, and a
memory-budget analysis that anchors the largest full solve the 8 GB box can
carry.

The current push is the unstructured path: a research note steers the
cell-centered finite-volume prototype, verified on a quad mesh, and a
Spalart-Allmaras research note scopes nuFor's first turbulence model. These
pages are the working record of that trajectory.
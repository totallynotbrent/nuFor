---
title: nuFor
---

nuFor is a CPU-first computational fluid dynamics code: finite-volume solvers
for the compressible Euler and Navier-Stokes equations, written as Fortran
numerical kernels behind a Rust application layer, with Python tooling for
verification. It runs structured cases in 1D, 2D, and 3D, from the terminal or
a web UI, with every numerical claim backed by a verification case.

## New here? Start with these

- [Getting started](getting-started.md) - install, run your first case, see a result
- [Architecture](architecture.md) - how the Rust/Fortran/Python layers fit together
- [Case format](formats/case-toml.md) - the `case.toml` file that drives every run

## Solver trajectory

**1D Euler** → **2D (Euler + viscous)** → **3D Euler**, then the unstructured
prototype and turbulence research. Follow the progression:

- [1D Euler](numerics/1d/1d-euler.md) - the finite-volume core: HLL flux, CFL, boundaries
- [2D Euler](numerics/2d/2d-euler.md) - HLLC, MUSCL reconstruction, the van Leer limiter
- [Viscous terms](numerics/2d/viscous.md) - the Navier-Stokes diffusive flux
- [3D solver](numerics/3d/3d-solver.md) - the structured 3D HLLC step
- [Unstructured FV](numerics/unstructured/unstructured-fv.md) - the research prototype

## Verification

Every milestone closes against an analytical or benchmark case:

- [Sod](numerics/1d/sod.md) and [Lax](numerics/1d/lax.md) - 1D shock tubes vs the exact Riemann solution
- [Oblique shock](numerics/2d/oblique-shock.md) - supersonic wedge vs θ-β-M theory
- [Shock reflection](numerics/2d/shock-reflection.md) - two-shock reflection off a wall
- [Channel flow](numerics/2d/channel-flow.md) - Poiseuille's parabolic profile
- [Supersonic cylinder](numerics/2d/supersonic-cylinder.md) - the detached bow shock
- [Regression suite](verification/regression-suite.md) - the automated gate on every change

## Reference

- [Equation of state](numerics/foundations/eos.md)
- [HLL flux](numerics/fluxes/flux-hll.md) · [HLL vs HLLC](numerics/fluxes/hll-vs-hllc.md)
- [State and grid](numerics/foundations/state-grid.md) · [Time step](numerics/foundations/time-step.md)
- [Restart](numerics/foundations/restart.md) · [HDF5 output](numerics/foundations/hdf5.md)
- [Threading](numerics/parallel/threading.md) · [SIMD](numerics/parallel/simd.md)

## Status

v1.2.0 (2026-09-09) is the current build: 1D/2D/3D structured Euler, 2D
viscous terms, CPU threading and AVX2 kernels, case-driven runs, and a web UI
with live field views. Open work includes turbulence implementation, a general
unstructured solver, thermochemistry, MPI, and broad experimental validation.
That remaining scope is tracked in the repository. See the
[full changelog](https://github.com/totallynotbrent/nuFor/blob/main/CHANGELOG.md).
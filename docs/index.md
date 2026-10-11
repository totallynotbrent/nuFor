---
title: nuFor
---

> **Work in progress.** nuFor is not complete. The validated core today:
> 1D/2D/3D Euler, 2D/3D viscous Navier-Stokes, Spalart-Allmaras RANS,
> masked bodies and cut cells, an axisymmetric solver, an equilibrium-air
> closure (CEA table + TGAS1 fits) on the hypersonic 2D paths, and an
> MPI-parallel Fortran march for the 2D equilibrium-air cases,
> exercised across a four-station atmospheric-entry sweep with a
> 1280x640 anchor run. Spalart-Allmaras 3D, two-temperature
> nonequilibrium, and unstructured meshes are the next rungs on the
> [roadmap](roadmap.md). Cases run from the CLI today; the web UI
> covers the fixed case set with field visualization, probes, and
> CSV / VTK / HDF5 export.

nuFor is a CPU-first computational fluid dynamics code: finite-volume solvers
for the compressible Euler and Navier-Stokes equations, written as Fortran
numerical kernels behind a Rust application layer, with Python tooling for
verification. It runs structured cases in 1D, 2D, and 3D, from the terminal or
a web UI, with every numerical claim backed by a verification case.

## New here? Start with these

- [Getting started](getting-started.md) - install, run your first case, see a result
- [Architecture](architecture.md) - how the Rust/Fortran/Python layers fit together
- [Case format](formats/case-toml.md) - the `case.toml` file that drives every run
- [Immersed bodies](formats/case-body.md) - the `[body]` section, axisymmetric runs, the boundary kinds

## Solver trajectory

**1D Euler** → **2D (Euler + viscous)** → **3D Euler**, then the unstructured
prototype and turbulence research. Follow the progression:

- [1D Euler](numerics/1d/1d-euler.md) - the finite-volume core: HLL flux, CFL, boundaries
- [2D Euler](numerics/2d/2d-euler.md) - HLLC, MUSCL reconstruction, the van Leer limiter
- [Viscous terms](numerics/2d/viscous.md) - the Navier-Stokes diffusive flux
- [3D solver](numerics/3d/3d-solver.md) - the structured 3D HLLC step
- [3D Navier-Stokes](numerics/3d/3d-viscous.md) - the viscous operator and no-slip walls in 3D
- [Unstructured FV](numerics/unstructured/unstructured-fv.md) - the research prototype

The full ladder, with what is explicitly out of scope, lives in the
[roadmap](roadmap.md).

## Research

- [Research results](research/results.md) - measured numbers behind the features: the station sweep, the trim anchor, MPI scaling, threading

## Verification

Every milestone closes against an analytical or benchmark case:

- [Sod](numerics/1d/sod.md) and [Lax](numerics/1d/lax.md) - 1D shock tubes vs the exact Riemann solution
- [Oblique shock](numerics/2d/oblique-shock.md) - supersonic wedge vs θ-β-M theory
- [Shock reflection](numerics/2d/shock-reflection.md) - two-shock reflection off a wall
- [Channel flow](numerics/2d/channel-flow.md) - Poiseuille's parabolic profile
- [Blasius flat plate](numerics/2d/blasius.md) - the laminar layer held steady, skin friction vs the exact correlation
- [Turbulent flat plate](numerics/2d/turbulent-plate.md) - the SA layer at correlation Reynolds numbers
- [Supersonic cylinder](numerics/2d/supersonic-cylinder.md) - the detached bow shock
- [Axisymmetric Euler](numerics/2d/axisymmetric.md) - bodies of revolution: the annular update and its geometric source
- [Axisymmetric sphere validation](verification/axisymmetric-sphere.md) - the bow shock standoff vs Ambrosio-Wortman
- [Aeroshell cross-check](aeroshell.md) - the capsule at peak q: CFD vs modified Newtonian
- [Regression suite](verification/regression-suite.md) - the automated gate on every change
- [Grid convergence study](verification/grid-convergence.md) - the turbulent plate refined 2x in every cell
- [3D channel validation](verification/viscous3d-channel.md) - the Poiseuille parabola held in a 3D slab
- [Square duct validation](verification/viscous3d-duct.md) - the four-wall series solution held steady

## Reference

- [Equation of state](numerics/foundations/eos.md) · [Equilibrium air](numerics/foundations/equilibrium-air.md)
- [HLL flux](numerics/fluxes/flux-hll.md) · [HLL vs HLLC](numerics/fluxes/hll-vs-hllc.md)
- [State and grid](numerics/foundations/state-grid.md) · [Time step](numerics/foundations/time-step.md)
- [Restart](numerics/foundations/restart.md) · [HDF5 output](numerics/foundations/hdf5.md)
- [Threading](numerics/parallel/threading.md) · [SIMD](numerics/parallel/simd.md) · [MPI](numerics/parallel/mpi.md)

## Status

The current build adds the MPI-parallel Fortran march to the 1.5.0 base:
1D/2D/3D structured Euler, 2D viscous terms on the equilibrium-air
closure, CPU threading and AVX2 kernels, case-driven runs, and a web UI
with live field views. Open work includes turbulence implementation, a
general unstructured solver, thermochemistry, and broad experimental
validation. That remaining scope is tracked in the
[full changelog](https://github.com/totallynotbrent/nuFor/blob/main/CHANGELOG.md).
---
title: nuFor
---

> **Work in progress.** nuFor is not complete. The validated core today:
> 1D/2D/3D Euler, 2D/3D viscous Navier-Stokes, Spalart-Allmaras RANS,
> masked bodies and cut cells, an axisymmetric solver, an equilibrium-air
> closure (CEA table + TGAS1 fits) on the hypersonic 2D paths, and an
> MPI-parallel Fortran march for the 2D equilibrium-air cases,
> exercised across a four-station atmospheric-entry sweep with a
> 1280x640 anchor run, and live two-temperature nonequilibrium
> (Park 2T) on the research path, plus the wall-resolved stagnation
> heating gate (Fay-Riddell, rung 6). Species continuity and
> unstructured meshes are the next rungs on the
> [[roadmap|roadmap]]. Cases run from the CLI today; the web UI
> covers the fixed case set with field visualization, probes, and
> CSV / VTK / HDF5 export.

nuFor is a CPU-first computational fluid dynamics code I've been building
for my hypersonic re-entry work: finite-volume solvers for the compressible
Euler and Navier-Stokes equations, written as Fortran numerical kernels
behind a Rust application layer, with Python tooling for verification. It
runs structured cases in 1D, 2D, and 3D, from the terminal or a web UI, with
every numerical claim backed by a verification case.

## New here? Start with these

- [[getting-started|Getting started]] - install, run your first case, see a result
- [[architecture|Architecture]] - how the Rust/Fortran/Python layers fit together
- [[formats/case-toml|Case format]] - the `case.toml` file that drives every run
- [[formats/case-body|Immersed bodies]] - the `[body]` section, axisymmetric runs, the boundary kinds

## Solver trajectory

**1D Euler** → **2D (Euler + viscous)** → **3D Euler**, then the unstructured
prototype and turbulence research. Follow the progression:

- [[numerics/1d/euler|1D Euler]] - the finite-volume core: HLL flux, CFL, boundaries
- [[numerics/2d/euler|2D Euler]] - HLLC, MUSCL reconstruction, the van Leer limiter
- [[numerics/2d/viscous|Viscous terms]] - the Navier-Stokes diffusive flux
- [[numerics/3d/solver|3D solver]] - the structured 3D HLLC step
- [[numerics/3d/viscous|3D Navier-Stokes]] - the viscous operator and no-slip walls in 3D
- [[numerics/unstructured/unstructured-fv|Unstructured FV]] - the research prototype

The full ladder, with what is explicitly out of scope, lives in the
[[roadmap|roadmap]].

## Research

- [[research/results|Research results]] - measured numbers behind the features: the station sweep, the trim anchor, MPI scaling, threading

## Verification

Every milestone closes against an analytical or benchmark case, tracked on the
[[verification/index|verification ladder]]:

- [[numerics/1d/sod|Sod]] and [[numerics/1d/lax|Lax]] - 1D shock tubes vs the exact Riemann solution
- [[numerics/2d/oblique-shock|Oblique shock]] - supersonic wedge vs θ-β-M theory
- [[numerics/2d/shock-reflection|Shock reflection]] - two-shock reflection off a wall
- [[numerics/2d/channel-flow|Channel flow]] - Poiseuille's parabolic profile
- [[numerics/2d/blasius|Blasius flat plate]] - the laminar layer held steady, skin friction vs the exact correlation
- [[numerics/2d/turbulent-plate|Turbulent flat plate]] - the SA layer at correlation Reynolds numbers
- [[numerics/2d/supersonic-cylinder|Supersonic cylinder]] - the detached bow shock
- [[numerics/2d/axisymmetric|Axisymmetric Euler]] - bodies of revolution: the annular update and its geometric source
- [[verification/axisymmetric-sphere|Axisymmetric sphere validation]] - the bow shock standoff vs Ambrosio-Wortman
- [[aeroshell|Aeroshell cross-check]] - the capsule at peak q: CFD vs modified Newtonian
- [[verification/regression-suite|Regression suite]] - the automated gate on every change
- [[verification/grid-convergence|Grid convergence study]] - the turbulent plate refined 2x in every cell
- [[verification/viscous3d-channel|3D channel validation]] - the Poiseuille parabola held in a 3D slab
- [[verification/viscous3d-duct|Square duct validation]] - the four-wall series solution held steady

## Reference

- [[numerics/foundations/eos|Equation of state]] · [[numerics/foundations/equilibrium-air|Equilibrium air]]
- [[numerics/fluxes/hll|HLL flux]] · [[numerics/fluxes/hll-vs-hllc|HLL vs HLLC]]
- [[numerics/foundations/state-grid|State and grid]] · [[numerics/foundations/time-step|Time step]]
- [[numerics/foundations/restart|Restart]] · [[numerics/foundations/hdf5|HDF5 output]]
- [[numerics/parallel/threading|Threading]] · [[numerics/parallel/simd|SIMD]] · [[numerics/parallel/mpi|MPI]]

## Status

The current build: 1D/2D/3D structured Euler, 2D/3D viscous terms,
Spalart-Allmaras RANS, masked bodies, an axisymmetric solver, the
equilibrium-air closure, live two-temperature nonequilibrium (Park 2T),
five-species finite-rate chemistry (Park 5), and an MPI-parallel
Fortran march. Open work: unstructured meshes, catalytic walls, and
broad experimental validation. The
[[roadmap|roadmap]] tracks it rung by rung.
---
title: nuFor — docs
description: nuFor, a CPU-first computational fluid dynamics framework.
---
# nuFor

nuFor is a portable, open-source, CPU-first computational fluid dynamics
framework: finite-volume solvers for the Euler and Navier-Stokes equations in
a Rust application layer, with Fortran numerical kernels and a Python analysis
stack. This site documents the architecture, numerics, and internal formats.

## Explore

- [[architecture|Architecture]] — system layout, tech baseline, build
- [[1d-euler|Numerics — 1D Euler milestone]] — current solver focus
- [[state-grid|Numerics — 1D state and uniform grid]] — state layout and mesh
- [[eos|Numerics — Ideal-gas equation of state]] — thermodynamics and pressure
- [[numerics/flux-hll|Numerics — HLL flux]] — the baseline approximate solver
- [[numerics/time-step|Numerics — CFL time-step control]] — the explicit step
- [[case-toml|Formats — case.toml]] — the reproducibility root of a run
- [[ffi-boundary|Research — FFI boundary]] — how Rust calls Fortran
- [[toml-config|Research — config crate choice]] — why TOML + `toml`
- [[research/flux-hll|Research — HLL vs HLLC]] — baseline flux decision
- [[research/cfl-time-step|Research — time-step rule]] — the CFL decision

## Status

Foundations are in place: mixed Rust/Fortran build, case schema and
validation, the FFI boundary, the 1D state and uniform grid geometry, the
ideal-gas equation of state, the HLL baseline flux, and CFL time-step control
are landed and tested. The next milestone is the 1D Euler solver itself; see
the roadmap in [[1d-euler|Numerics]] and the broader plan in the repository.
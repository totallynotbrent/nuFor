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
- [[case-toml|Formats — case.toml]] — the reproducibility root of a run
- [[ffi-boundary|Research — FFI boundary]] — how Rust calls Fortran
- [[toml-config|Research — config crate choice]] — why TOML + `toml`

## Status

Foundations are in place: mixed Rust/Fortran build, case schema and
validation, and the FFI boundary are landed and tested. The next milestones
build the 1D and then 2D Euler solvers; see the roadmap in
[[1d-euler|Numerics]] and the broader plan in the repository.
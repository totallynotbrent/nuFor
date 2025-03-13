# nuFor

A CPU-first computational fluid dynamics research code. nuFor couples modern
Fortran numerical kernels with a Rust application layer (CLI, web UI, runtime)
and Python tooling for analysis, reference solutions, and verification.

The project is being developed openly over roughly five years, one verified step
at a time: foundations and 1D Euler first, then 2D, viscosity, parallel
performance, and broader validation. This repository is the executable record of
that program. Every numerical claim is meant to be backed by verification
evidence, not by a successful compile.

## Status

Experimental. The mixed Rust/Fortran build system and the FFI boundary are in
place (PLAN step 2): a Cargo workspace drives CMake via `build.rs`, the
Fortran kernel library compiles and links, and FFI integration tests pass.
The boundary decision is recorded in `docs/research/ffi-boundary.md`. Numerical
solvers come next. See `plans/PLAN.md` for the ordered roadmap and current step.

## Repository layout

```
nuFor/
├── README.md            this file
├── LICENSE              GPL-3.0-only
├── Cargo.toml           Rust workspace
├── CMakeLists.txt       (future top-level native orchestration)
├── CONTRIBUTING.md      how to contribute
├── CODE_OF_CONDUCT.md   contributor expectations
├── SECURITY.md          security reporting policy
├── docs/                architecture, numerics, physics, formats, research,
│                        verification, operations
├── plans/               step-by-step roadmap and execution log
├── crates/              Rust packages (nufor-core, then nufor-cli/runtime/...)
├── fortran/             Fortran numerical kernels (CMake build)
├── python/              analysis, reference, visualization, verification
├── cases/               tutorial, verification, benchmark cases
├── tests/               unit, integration, regression, cross-platform
├── tools/               project tooling
└── scripts/             helper scripts
```

The full engineering constitution (spec), research notebook index, and roadmap
live in `nuFor_README.md`. It is the authority on architecture, numerical
standards, and milestone gates.

## Toolchain

The dev reference environment uses:

| Component | Version |
|-----------|---------|
| gfortran  | 14.3 (GCC 14) |
| Rust      | 1.98 (cargo/rustc) |
| CMake     | 3.31 |
| HDF5      | 1.14 |

Builds target Linux; Windows support is tracked in CI.

## Building

`cargo test --workspace` builds and tests the full stack: `build.rs` drives
CMake to compile the Fortran kernels, links the static archive, and runs FFI
integration tests. To build the Fortran kernels standalone (see also
`docs/operations/README.md`):

```
cmake -S fortran -B build/fortran -DCMAKE_BUILD_TYPE=Release
cmake --build build/fortran --config Release
```

## Roadmap

- Year 1: foundations, 1D Euler, verification, basic architecture, first web/CLI shell.
- Year 2: 2D Euler, better reconstruction/fluxes, performance foundation, 2D viscous beginnings.
- Year 3: 3D, unstructured mesh research, turbulence foundations, stronger visualization.
- Year 4: RANS/thermo/propulsion workflows, MPI research, mature data formats, benchmarks.
- Year 5: integration, broad validation, documentation, release quality, publication readiness.

## License

GPL-3.0-only. See `LICENSE`.
# nuFor

A CPU-first computational fluid dynamics research code. nuFor couples modern
Fortran numerical kernels with a Rust application layer (CLI, web UI, runtime)
and Python tooling for analysis, reference solutions, and verification.

The project is developed one verified step at a time: foundations and the 1D
Euler solver first, then 2D, viscosity, parallel performance, and broader
validation. This repository is the executable record of that program. Every
numerical claim is meant to be backed by verification evidence, not by a
successful compile.

## Status

**Stable public-ready build, v1.2.0.** The core solver path is verified
end-to-end: mixed Rust/Fortran build, 1D/2D/3D Euler, HLL/HLLC fluxes, MUSCL
reconstruction, 2D viscous terms with channel-flow verification, CPU
threading and AVX2 kernels, restart files, case-file driven runs, and a web
UI with live field and 3D-slice views. 34+ test suites gate every commit.

Not yet: turbulence implementation, unstructured meshes beyond the research
prototype, thermochemistry/reacting flow, MPI, and the broad experimental
validation database a mature solver carries.

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

## Running the web UI

```
cargo run --release -- serve 8060
```

The server binds `0.0.0.0:8060`, so open `http://<host>:8060` from the same
machine, over LAN, or via Tailscale. Use the release build — the debug binary
is an order of magnitude slower at the flow solves behind the UI. The page is
a single CFD workspace: a central flow viewport (2D blast animation or 3D
slice planes, with a toggleable mesh overlay and colorbar), a left
data/display sidebar, and a diagnostics dock (1D monitor, line probes,
comparison overlay, run history). The Run button executes your configured case
(2D/3D), saves `results/<name>.vtk`, and reports the convergence summary.

## Running a case

One `case.toml` drives the solver from the terminal across all three dimensions:

```
nufor init case.toml              # write a runnable template
nufor run case.toml               # solve it (euler_1d / euler_2d / euler_3d)
```

`physics.equations` picks the dimension, `[mesh]` the cell count and domain
(2D/3D add `ny/nz`), and `initial_condition` is `uniform`, `two_state` (1D), or
`blast` (2D/3D). Results land next to the case as csv/vtk (and a PNG for 2D).
To bring your own mesh, set `mesh.source = "file"` with a `path` to a list of
cell-center coordinates (must be uniformly spaced). Schema and validation in
`docs/formats/case-toml.md`.

To keep the web UI running unattended, run the server under a systemd user
unit so it stays up and auto-restarts on rebuild.

## Roadmap

Done through v1.2.0 (see CHANGELOG):
- 1D/2D/3D structured Euler solvers with case-file–driven runs
- 2D viscous terms, MUSCL + HLL/HLLC, positivity, CPU threading + AVX2
- Restart files, HDF5/VTK/CSV output, line probes, comparison dashboard
- Web UI: flow viewport, mesh overlay, 3D slice view, case runs from the page

What's left, roughly in order:
- **Viscous maturity** — laminar validation against experiment, transport-model
  docs, verified wall boundary behavior (3D is Euler-only today)
- **Unstructured FV** — the prototype on quads exists; a general unstructured
  solver with real mesh import (Gmsh/VTK) does not
- **Turbulence** — Spalart–Allmaras is at the research-note stage; needs
  implementation, verification, wall treatment, and benchmarks
- **Thermodynamic property expansion + reacting flow** — species, chemistry
  source terms, high-temperature gas effects (the propulsion scope)
- **HPC maturity** — MPI/domain decomposition, scaling studies, reproducible
  benchmark harness (Intel + AMD), NUMA awareness
- **Mesh ecosystem** — geometry import (STL/STEP), mesh-quality diagnostics,
  boundary tagging
- **Validation database** — broad experimental/benchmark comparisons beyond
  the current verification cases
- **Year-5 polish** — cross-platform audit, license/security audit, first
  stable public release, published benchmarks

## License

GPL-3.0-only. See `LICENSE`.
---
title: Architecture
---

nuFor separates a Rust application layer from a Fortran numerical layer across
a narrow C ABI, with Python reserved for analysis and reference calculations
(never production hot loops).

```mermaid
graph TD
    UI[Browser / web UI] -->|HTTP + WebSocket| RUST[Rust app layer]
    RUST[CLI, web API, config, case/job mgmt, logging]
    RUST -->|stable FFI boundary| FTN[Fortran numerical layer]
    FTN[EOS, fluxes, reconstruction, Riemann solvers, viscous, kernels]
    FTN --> OMP[OpenMP]
    FTN --> CPU[Intel / AMD CPU]
```

## Language split

- **Rust**: application framework, CLI, web/API, orchestration, integration.
- **Fortran**: numerical and physics kernels and dense array operations.
- **Python**: analysis, reference solutions, plotting, verification,
  automation.

## Interop

Rust and Fortran cross a deliberate C-compatible ABI: contiguous arrays with
explicit lengths, structured error codes, no Fortran derived-type leakage.
See [[numerics/foundations/memory-layout|the memory layout]] for the full design and its rules.

## Build

Cargo owns the Rust workspace; a `build.rs` drives CMake to compile the
Fortran kernels into a static archive and link them in. The tree is
`crates/` (Rust), `fortran/` (Fortran), `python/` (analysis), `cases/`
(runnable cases), `docs/` (this site).

## Areas

- [[numerics/1d/euler|Numerics - 1D Euler]]
- [[formats/case-toml|Formats - case.toml]]

## See also

- [[numerics/tools/web-ui|Web UI]]
- [[formats/case-toml|case.toml format]]
- [[numerics/3d/state-grid|3D foundations]]
- [[performance/optimization-ledger|Optimization ledger]]

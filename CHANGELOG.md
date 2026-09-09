# Changelog

All notable changes to nuFor are tracked here. This project follows a
graduated, verification-first release cadence.

## [1.1.0] - 2026-09-09

A feature release focused on the way nuFor is looked at: the web UI is now a
single CFD-style workspace rather than a page of tabs, and it can show 3D
fields.

### Web UI

- Rebuilt as one workspace like ParaView/Tecplot: a central flow viewport
  with a data/display sidebar on the left and a diagnostics dock below. No
  per-feature tabs.
- Viewport shows the solved flow field with an optional toggleable cartesian
  mesh overlay, a pinned colorbar, and live play/scrub over the blast time.
- 3D mode: slice planes through a resolved 3D sphere-blast, with the
  coordinate axis and slice position selectable.
- Left sidebar holds the simulation controls (case, cells, t, gamma, cfl,
  run) and display controls (dimension, field, layers, mesh density,
  playback). The dock holds the 1D monitor, line probes, the comparison
  overlay, and run history.
- The visualization scales to the window; it is served on `0.0.0.0` so it is
  reachable over LAN or Tailscale.

### Running the web UI

```
cargo run --release -- serve 8060
```

then open `http://<host>:8060`. Use the release binary — the flow solves are
orders of magnitude faster than the debug build.

## [1.0.0] - 2026-09-05

The first stable public release. The solver covers the 1D, 2D, and 3D Euler
equations (and the 2D Navier-Stokes terms) in a Rust application shell with a
Fortran-open, CPU-first design, shipped with a live Quartz documentation site,
a web dashboard, and a verification suite that reproduces the milestone
results.

### Numerics

- 1D finite-volume Euler solver: HLL and HLLC fluxes, CFL time-step control,
  MUSCL reconstruction with the van Leer limiter, verified against the exact
  Riemann solution (Sod, Lax, and the four-rung ladder).
- 2D Euler solver with HLLC flux, oblique-shock and shock-reflection cases
  verified against shock relations, positivity-preserving update.
- 2D viscous Navier-Stokes terms, Poiseuille-channel validation.
- 3D structured solver with a 3D HLLC step, verified against the Sod tube,
  and a documented memory budget against the 8 GB target.
- Unstructured cell-centered prototype, verified on a quad mesh; Spalart-Allmaras
  research note scopes the first turbulence model.
- Supersonic-cylinder immersed body verified against a detached bow shock
  reproducing the normal-shock density ratio.
- Runtime-dispatched AVX2 kernels in the conservation update, with a scalar
  fallback and full numerical-equivalence tests.
- Restart format hardened with an explicit schema version, strict size and
  physical-sanity validation, and legacy reads.

### Interfaces

- Web UI with live 1D plots, 2D field images, line probes, and a comparison
  dashboard overlaying fields and resolutions.
- CSV, VTK, and HDF5 output plus a bit-exact binary restart.
- CLI: `run`, `inspect`, `export`, `benchmark`, `serve`, `version`.
- Quartz documentation site with a page graph, mermaid, and LaTeX, published
  on GitHub Pages.

### Verification

Cross-platform CI (Linux and Windows) runs the Rust build, the Fortran-linked
FFI integration tests, and a skeleton/LF gate on every change. The full
workspace suite is reproducible from a clean checkout.
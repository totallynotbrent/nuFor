# Changelog

## [1.4.0] - 2026-09-21

Turbulence: the solver now carries the Spalart-Allmaras model. A case with
`equations = "rans_2d_sa"` solves the 2D Navier-Stokes equations coupled to
the transported SA variable, with the eddy viscosity feeding the momentum
and heat fluxes. The model follows the 1994 baseline (no trip terms) with
the standard guards, closes at `nu_tilde = 0` exactly, and reproduces the
laminar solver bit for bit in that limit.

No-slip walls now shear the flow. The viscous pass carries a
quadratic-consistent one-sided wall flux, so a boundary layer actually
develops against a wall; the Poiseuille hold balances to machine precision
on the discrete stencil. `nufor run` prints skin friction along solid walls
after an SA case, and `cases/flat-plate/case.toml` ships as the runnable
plate example.

Also new: `[boundaries]` gained `top`/`bottom` sides for the 2D solvers, and
`[physics]` gained `mu`, `pr`, and a `[physics.turbulence]` block.

## [1.3.0] - 2026-09-14

Mesh diagnostics: `nufor mesh-check case.msh` validates an imported
unstructured mesh before solving. It reports cell/face counts, area
min/mean/max and the stretch ratio, a per-cell closure residual, flipped
interior faces, and negative or degenerate volumes, then gives a verdict.
The same summary is available at `GET /api/mesh-stats?path=...` for tooling.

`Ugrid` now stores cell centroids, which feeds the orientation and closure
checks.

## [1.2.0] - 2026-09-09

`nufor run <case.toml>` now dispatches on `physics.equations` and solves 1D
(`euler_1d`), 2D (`euler_2d`), or 3D (`euler_3d`) from a single case file.
Results (csv/vtk, plus a PNG for 2D) are written next to the case and a
summary is printed.

- `[mesh]` adds the 2D/3D edges (`ny/y0/y1`, `nz/z0/z1`). `source = "file"`
  plus `path` loads cell-center coordinates from a file. Coordinates must be
  uniformly spaced; non-uniform input is refused with a clear message.
- New `blast` initial condition, an over-pressured fireball in ambient
  surroundings, for 2D/3D. `uniform` and `two_state` stay for 1D. `nufor init`
  writes a runnable template. Validation rejects a non-positive blast radius.
- `write_vtk2d` and `write_vtk3d` emit ParaView-readable `STRUCTURED_POINTS`
  blocks (density/pressure/energy).
- The web Run button drives `/api/run-case` with the configured dimension,
  cell count, final time, and CFL, saves `results/<name>.vtk`, and shows the
  convergence summary (steps, time, wall clock) in the Monitor dock.
- The UI header version now reads from the crate instead of a hardcoded string.

## [1.1.0] - 2026-09-09

### Web UI

- Rebuilt as one workspace like ParaView/Tecplot: a central flow viewport, a
  data/display sidebar on the left, a diagnostics dock below. No per-feature
  tabs.
- Viewport shows the solved field with a toggleable cartesian mesh overlay, a
  pinned colorbar, and play/scrub over the blast time.
- 3D mode: slice planes through a solved 3D sphere-blast, with the axis and
  slice position selectable.
- The sidebar holds simulation controls (case, cells, t, gamma, cfl, run) and
  display controls (dimension, field, layers, mesh density, playback). The
  dock holds the 1D monitor, line probes, the comparison overlay, and run
  history.
- Served on `0.0.0.0` so it is reachable over LAN or Tailscale.

### Running the web UI

```
cargo run --release -- serve 8060
```

then open `http://<host>:8060`. Use the release binary; the flow solves are
orders of magnitude faster than the debug build.

## [1.0.0] - 2026-09-05

First stable public release. The solver covers the 1D, 2D, and 3D Euler
equations and the 2D Navier-Stokes terms in a Rust application shell with a
CPU-first design. It ships with a Quartz documentation site, a web dashboard,
and a verification suite.

### Numerics

- 1D finite-volume Euler solver: HLL and HLLC fluxes, CFL time-step control,
  MUSCL reconstruction with the van Leer limiter, verified against the exact
  Riemann solution (Sod, Lax, and the four-rung ladder).
- 2D Euler solver with HLLC flux. Oblique-shock and shock-reflection cases
  verified against shock relations. Positivity-preserving update.
- 2D viscous Navier-Stokes terms, Poiseuille-channel validation.
- 3D structured solver with a 3D HLLC step, verified against the Sod tube, and
  a documented memory budget against the 8 GB target.
- Unstructured cell-centered prototype, verified on a quad mesh. A
  Spalart-Allmaras research note scopes the first turbulence model.
- Supersonic-cylinder immersed body verified against a detached bow shock
  reproducing the normal-shock density ratio.
- Runtime-dispatched AVX2 kernels in the conservation update, with a scalar
  fallback and numerical-equivalence tests.
- Restart format with an explicit schema version, size and physical-sanity
  validation, and legacy reads.

### Interfaces

- Web UI with live 1D plots, 2D field images, line probes, and a comparison
  dashboard.
- CSV, VTK, and HDF5 output plus a bit-exact binary restart.
- CLI: `run`, `inspect`, `export`, `benchmark`, `serve`, `version`.
- Quartz documentation site with a page graph, mermaid, and LaTeX, published
  on GitHub Pages.

### Verification

Cross-platform CI (Linux and Windows) runs the Rust build, the Fortran-linked
FFI integration tests, and a skeleton/LF gate on every change. The full
workspace suite is reproducible from a clean checkout.
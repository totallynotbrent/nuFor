---
title: Web UI skeleton
---

A dependency-free HTTP server in the `nufor` CLI (`nufor serve`) exposes the
whole solvable surface as a small data API plus an app shell page. The front
end is intentionally plain and functional: an external design tool restyles it
against the same endpoints, so the API is the contract that stays stable.

## Running

    nufor serve [port]        # default port 8060, bound to 0.0.0.0

The server listens on all interfaces so the dashboard is reachable over the
LAN or Tailscale, not just on the box itself.

## Pages

The `/` shell has one panel per feature:

- **Solve** — pick the field (density / velocity / pressure / Mach / momentum /
  energy), run, and watch the curve; optional exact-solution overlay.
- **Case** — edit the case parameters used by the next run (initial condition,
  cells, end time, gamma, CFL, boundary).
- **Verify** — run against the exact 1D Riemann solution and read the L1
  density error.
- **Output** — download the current snapshot as CSV, VTK, or HDF5.
- **Benchmark** — run a single-core throughput sweep.
- **History** — every run in the session with its mesh, outcome, and reason.
- **Gallery** — the validated case set rendered beside its mesh: the 2D
  blast, the laminar Blasius plate, and the turbulent SA plate. The mesh
  pane shows the actual grid (uniform for the blast, wall-clustered for the
  plates) so the boundary-layer spacing is visible. The turbulent plate
  solves once in the background at server start and reports its progress
  until it lands.
- **2D field** — a 2D blast wave solved and rendered to a colour-mapped image,
  switchable between density, mach, and pressure. Playback animates from
  frames baked once in the background; the run itself marches in a
  background thread and reports live progress (percent, step count, and
  time) under the run button instead of blocking the page.
- **Validation** — the measured skin friction along the plate against the
  published reference: the exact Blasius correlation for the laminar plate,
  the Schlichting power law for the turbulent one.
- **Mesh view** — upload a mesh (2D or 3D rectilinear VTK / coordinate file,
  or a 2D gmsh `.msh`), see it on the left of the main screen with the
  result beside it on the right. Rectilinear meshes can be solved directly
  from the browser: a blast marches on the imported grid in the background
  and the finished field can be switched between density, mach, and
  pressure; 3D meshes add axis switching (xy / xz / zy) and a slice slider.
  3D meshes render as a drag-to-spin wireframe for inspection. gmsh meshes
  are wireframe-view only for now — solving them needs the unstructured
  path, which is CLI-driven today.

## Data API

All routes answer on `0.0.0.0:<port>`:

- `GET /api/result` (alias `/api/snapshot`) — the current result envelope:
  steps, residual, termination reason, the snapshot (centers, rho, m, e, u, p,
  mach), and the exact solution with its L1 density error.
- `GET /api/config` — the current case configuration as JSON.
- `GET /api/run?kind=sod|lax&n=..&t=..&gamma=..&cfl=..&bc=transmissive|reflective`
  — reconfigure and run, returning the same envelope; the run is recorded.
- `GET /api/exact` — just the exact arrays (rho, u, p) plus the L1 error.
- `GET /api/export?format=csv|vtk|h5` — the snapshot written to that format and
  returned as the file bytes.
- `GET /api/benchmark?steps=..` — a JSON throughput table across mesh sizes.
- `GET /api/history` — every run in the session as JSON.
- `GET /api/image?n=..&field=rho|mach|p` — a 2D blast wave solved at `n` cells
  and returned as a PNG of the chosen scalar field (`image/png`).
- `GET /api/run-case?dim=2d|3d&n=..&t=..&cfl=..&name=..` — start a background
  blast march; answers `202` immediately. Poll `GET /api/run-status` for
  percent, steps, and time; the finished run writes `results/<name>.vtk`.
- `GET /api/frames?t=..` — one playback frame (PNG) from the baked blast
  animation; `GET /api/frames/times` lists the available times.
- `GET /api/mesh-faces?kind=blast|laminar|turbulent` — the face coordinates of
  the active gallery grid, so the mesh pane draws the real cell spacing.
- `GET /api/plates/{laminar|turbulent}/image|validation|profile` — the plate
  gallery: a rendered field, the cf(x) table with its reference correlation,
  and the wall-normal u(y) profile. The turbulent variants answer `503`
  with a plain-text status until the background solve finishes;
  `GET /api/plates/progress` reports how far along it is.
- `POST /api/upload?name=..` — upload a mesh file (64 MB cap, 8M cell cap
  for the browser solve path); it lands in `uploads/` and gets parsed.
- `GET /api/meshview/wire` — the wireframe of the uploaded mesh: 2D face
  lines, a sampled 3D grid box for the spin view, or gmsh nodes + edges.
- `GET /api/meshview/solve?t=..&cfl=..` — start the background march on the
  uploaded rectilinear grid; poll `GET /api/meshview/status` for progress.
- `GET /api/meshview/slice?field=..&axis=..&u=..` — the solved field as a
  PNG: the full 2D field, or the plane at `u` along `axis` for 3D.

The rule going forward: whenever a new feature lands (2D HLLC, wedge cases,
viscous terms), it gets an API route and a panel in the same commit, so the
skeleton stays a faithful mirror of the solver.

## See also

- [[architecture|Architecture]]
- [[numerics/2d/2d||2D foundations]]
- [[numerics/1d/1d-euler||1D Euler]]
- [[numerics/tools/line-probes||Line probes]]
- [[numerics/tools/comparison-dashboard||Comparison dashboard]]

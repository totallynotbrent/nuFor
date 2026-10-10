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

The `/` shell is one solver workspace in the shape of a commercial CFD
workbench: a header with the case file actions, a left setup sidebar, and
one large graphics window that serves both the mesh and the result.

- **Header** — **Load case** opens a picker over the `cases/` directory
  (the list is the directory; nothing is baked into the page), **Save**
  writes the editor's changes back to the file, **New** saves the
  current setup under a fresh name, and **Run** marches the loaded case.
  The current case name and run status sit beside them, with a progress
  bar that fills as the march runs.
- **Setup** — the loaded case as a tree: Mesh, Body, Physics, Initial
  condition, Boundaries, Numerics, Time control, Output, Metadata, each
  node a set of property inputs (enum fields render as selects). Edits
  write into the case file on Save; "Edit raw TOML" flips to the full
  text editor.
- **Graphics window** — two modes on the window itself, Mesh and
  Results. **Mesh** renders the loaded case's grid like a mesher: cell
  edges, the body filled with its outline, each domain boundary
  colored by its condition with a legend naming them, a cell count,
  and wheel zoom / drag pan. **Results** shows the field (density /
  mach / pressure) at true aspect with coordinate axes, a scrub bar
  over the run's frames, and a mesh overlay toggle. Frames arrive
  while the run is still marching, so the scrub bar and the field can
  be watched live; a case with no frames yet shows a placeholder
  message instead of a broken image. For dimensionally small runs the
  displayed times round to 0.000 at three decimals.
- **Results panel** (sidebar) — after a run with a body: the surface Cp
  stations and the integrated axial-force coefficient, normalized by
  the case's own freestream and reference area.
- **Import mesh** (collapsed section) — upload or server-path load of
  a mesh (rectilinear VTK / coordinates, gmsh 2d); the uploaded grid
  renders in the same mesh viewer, and rectilinear grids can run a
  blast march on it.

## Data API

All routes answer on `0.0.0.0:<port>`:

- `GET /api/cases` — the case directory as a JSON list.
- `GET /api/cases/load?name=..` — one case's `case.toml` as text.
- `POST /api/cases/save?name=..` — write the editor's text back to the
  case file (validated before anything is written).
- `POST /api/cases/create?name=..` — save the editor's text as a new case.
- `POST /api/cases/delete?name=..` — remove a case directory.
- `GET /api/cases/run?name=..` — start the background march for that
  case (1D, 2D, or axisymmetric, dispatched by the case's equations);
  answers `202` immediately. Poll `GET /api/cases/run-status` for
  percent, steps, and time.
- `GET /api/cases/mesh?name=..` — the loaded case's grid faces and
  extents plus the body outline sampled from its own sdf, so the mesh
  pane shows the exact grid the march will run on.
- `GET /api/cases/probe?field=..&x0=..&y0=..&x1=..&y1=..&samples=..` —
  a line sample through the finished run's final field.
- `GET /api/cases/surface` — after a run with a body: the surface Cp
  stations, the body silhouette, and the integrated C_A (normalized by
  the case's own freestream and reference area).
- `GET /api/frames?t=..&field=rho|mach|p` — one playback frame (PNG)
  from the animation, per field; frames are served as they are
  rendered, during the run as well as after it.
  `GET /api/frames/times` lists the frame times, growing live.
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
- `GET /api/frames?t=..&field=rho|mach|p` — one playback frame (PNG)
  from the animation, per field; frames are served as they are
  rendered, during the run as well as after it.
  `GET /api/frames/times` lists the frame times, growing live.
- `GET /api/cases` — the case directory as a JSON list.
- `GET /api/cases/load?name=..` — one case's `case.toml` as text.
- `POST /api/cases/save?name=..` — write the editor's text back to the
  case file (validated before anything is written).
- `POST /api/cases/create?name=..` — save the editor's text as a new case.
- `POST /api/cases/delete?name=..` — remove a case directory.
- `GET /api/cases/run?name=..` — start the background march for that
  case; answers `202` immediately. Poll `GET /api/cases/run-status` for
  percent, steps, and time.
- `GET /api/cases/surface` — after a run with a body: the surface Cp
  stations, the body silhouette, and the integrated C_A (normalized by
  the case's own freestream and reference area).
- `GET /api/cases/mesh?name=..` — the loaded case's grid faces and
  extents plus the body outline sampled from its own sdf, so the mesh
  pane shows the exact grid the march will run on.
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

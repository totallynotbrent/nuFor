# Web UI skeleton

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
- **2D field** — a 2D blast wave solved and rendered to a colour-mapped image,
  switchable between density, mach, and pressure.

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

The rule going forward: whenever a new feature lands (2D HLLC, wedge cases,
viscous terms), it gets an API route and a panel in the same commit, so the
skeleton stays a faithful mirror of the solver.

## See also

- [[architecture|Architecture]]
- [[numerics/2d|2D foundations]]
- [[numerics/1d-euler|1D Euler]]

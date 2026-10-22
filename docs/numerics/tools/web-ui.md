---
title: Case console
---

The `nufor watch` command opens the live case console: one page that follows
a case from launch to finish and leaves a report behind. There is no editor
and no upload path, cases are written as TOML and run from the CLI, and the
console is a read-only window onto the march.

## Running

    nufor watch <case.toml> [port]   # standalone viewer, default port 8061
    nufor run <case.toml>            # boots the console automatically

Every `nufor run` spawns the console on a background thread before the
first step and prints the link up front:

    viewer: http://127.0.0.1:8062/  (live case console)

so the page is already open while the march populates it. The port scanner
takes the first free port in 8062..8070, several concurrent runs each get
their own console.

## Layout

The page is a three-column CAE console, full height, in the house design
(dark canvas, amber accent, Inter body, JetBrains Mono data, Newsreader
display).

- **Header**  -  case name, the equations line (1D Euler, 2D axisymmetric,
  RANS-SA, ...), the state pill (running / finished / failed), and a live
  progress box while a run is marching: a fill bar with the percentage of
  final_time covered and an ETA from the recent pace. The readout under
  it refreshes every step and the field image re-fetches on the poll
  cadence while the run is live, so the picture moves with the march.
- **Left rail**  -  the live run readout (state, step, physical time, dt,
  wall clock, sample count, march residual) plus the case card: grid,
  domain extents, body presence.
- **Center viewport**  -  tabbed figures, all sharing one physical fit so
  the mesh and every field sit in the same box at the same scale:
  - **field**  -  density, mach, pressure, temperature, and schlieren,
    each a button on the field bar. The first four render through the
    same renderer the results pipeline uses, body silhouette included;
    the color window comes from the data itself (2nd/98th percentile
    over live cells) and the colorbar shows it. Schlieren is the
    classic |grad rho| magnitude, log-scaled and grayscale: dark
    freestream, bright shock bands, the way experimental schlieren
    reads.
  - **mesh**  -  cell edges and the body outline at true aspect.
  - **forces**  -  C_D and C_L (or the density residual for cases without
    a body) as a history over the march.
  - **log**  -  the run log as real text, selectable and copyable,
    scrolled like any text pane instead of painted into the canvas.
- **Right rail**  -  the forces block (C_D, C_L, L/D, peak wall heat flux),
  the peaks block (p_max, T_max, rho_max, rho_min), the hover probe (cell
  indices, physical coordinates, and the local density, pressure,
  temperature, and mach under the cursor), the rho/p/dt history plot,
  and the force history (C_D, C_L, q_max over the march).
- **Left rail**  -  the run readout plus the full case card, from grid
  and domain through the closure, freestream Mach and state, wall
  temperature, flux and reconstruction order, cfl, and the stop
  conditions.

## Where the numbers come from

The march writes a heartbeat (`progress.json`, next to the case) every 50
steps, step 1 included. On those emit steps it also integrates the surface
force coefficients over the body band, the peak wall heat flux when the
case is viscous with a wall temperature, and the density residual. The
console polls that heartbeat once a second, so the readout, plots, and
force panels fill while the run is still marching.

The finished run's fields come from the written vtk: the console recovers
the state through the same reader the results pipeline uses and renders on
request, so the figures always match the case's own eos and model.

## The report

When a run finishes, `nufor run` also writes a standalone report next to
the case:

    report: /tmp/mycase/aeroshell-capsule-trim-report.html

One self-contained file: the console page with the progress history, run
log, mesh geometry, and every field render (density, mach, pressure,
temperature, schlieren) inlined as data urls. It opens in any browser with
no server, no network, and no nufor install, which makes it the artifact
to attach to a submission or print to pdf.

## Data API

`nufor watch` answers a small read-only API alongside the page:

- `GET /`  -  the console page.
- `GET /api/watch`  -  the heartbeat: kind, state, step, t, dt, samples
  (step, t, dt, rho/p extrema, wall, cd, cl, q_max, residual).
- `GET /api/mesh`  -  grid edges and extents plus the body outline.
- `GET /api/field?name=rho|mach|p|t|schlieren`  -  the finished field as a
  PNG, ranges taken from the data (404 before the first finished run).
- `GET /api/range?name=..`  -  the color window the colorbar reads.
- `GET /api/cell?i=&j=`  -  the finished field's local values (density,
  pressure, temperature, mach) at one cell, for the hover probe.
- `GET /api/case_info`  -  the case card: closure, freestream state and
  Mach, wall temperature, numerics, and stop conditions.
- `GET /api/log`  -  the run log text.

## See also

- [[architecture|Architecture]]
- [[numerics/2d/state-grid|2D foundations]]
- [[numerics/tools/line-probes|Line probes]]
- [[numerics/tools/comparison-dashboard|Comparison dashboard]]

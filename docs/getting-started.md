---
title: Getting started
---

This walks you from a clean checkout to a solved case in a few minutes.

## Prerequisites

- **Rust** (1.74+) and **cargo**
- **gfortran** (GCC 14) and **CMake** to build the Fortran kernels
- Optional: **HDF5** for the self-describing binary output

## Build and test

```bash
git clone https://github.com/totallynotbrent/nuFor.git
cd nuFor
cargo test --workspace
```

`cargo test` drives CMake to build the Fortran kernels, links them, and runs
the verification suite. A clean run means all the shock tubes, wedge, channel
flow, and regression cases pass on your machine.

## Run a case

Every run is defined by a single `case.toml`. Generate one, edit it, and run it:

```bash
cargo run --release --bin nufor -- init case.toml     # write a template
cargo run --release --bin nufor -- run  case.toml     # solve it
```

The template defaults to a 1D Sod shock tube. Change `physics.equations` to
`euler_2d` or `euler_3d`, set the mesh dimensions, and the same `run` command
solves the higher-dimensional blast case. Results (CSV, VTK, and a PNG for 2D)
are written next to the case file.

Every case run also writes two small logs next to the `case.toml`:
`run-log.txt` is a human-readable summary (grid, timestep count, wall clock,
freestream state, peak pressure and temperature), and `run-log.csv` appends
the same numbers as one row per run so a report or spreadsheet can diff
runs head-to-head.

| What you want | Key fields to change |
|---|---|
| A 2D/3D blast | `equations = "euler_2d"` or `"euler_3d"`, add `mesh.ny`/`nz`, `initial_condition.type = "blast"` |
| A custom mesh | `mesh.source = "file"` + `mesh.path` to a mesh file |
| A wall-clustered mesh | `mesh.source = "clustered"` + `first_cell`/`growth` (`cluster = "wall"`, `"channel"`, or `"top"`) |
| A boundary-layer inflow | `boundaries.left = "profile_inflow"` + `[boundaries.inflow_profile]` (a built-in Blasius layer or a `y u v` table) |
| A different end time | `time.final_time` |
| More resolution | `mesh.nx` (and `ny`/`nz`) |

The full schema is documented in [[formats/case-toml|the case format reference]].

Equilibrium-air 2D Euler cases can also run multi-rank. Build once, then
launch the same binary under `mpirun` (OpenMPI required):

```bash
mpirun -np 4 nufor run case.toml
```

np=1 takes the same code path, so results are identical either way; the
measured scaling and the parity gates are on [[numerics/parallel/mpi|the MPI page]].

## Watch it run

The console boots with every run. Launch a case and the viewer link prints
before the first step:

```bash
nufor run case.toml
# viewer: http://127.0.0.1:8062/  (live case console)
```

Open the link and the case console fills in live: the run readout (step,
time, dt, residual), the force coefficients as they converge, and, once
the march finishes, the field figures (density, mach, pressure,
temperature, schlieren) and the mesh. To open the console on an
already-finished case without rerunning:

```bash
nufor watch case.toml 8061
```

The layout and the data API behind it are on [[numerics/tools/web-ui|the
case console page]]. When the run finishes, a standalone HTML report lands
next to the case, one file with everything inlined that opens in any
browser with no server.

The console reads the case's files while the march is still writing
them, so it can catch a half-written field. A bad read is refused and
the run carries on, because a viewer glitch must never take a march
down with it.

## Bring your own mesh

For 2D/3D, set `mesh.source = "file"` with `mesh.path` pointing at a
**structured** mesh file - an ascii VTK `RECTILINEAR_GRID`/`STRUCTURED_POINTS`
block, or a plain `coordinates` file (one face coordinate per line, a `#` or
blank line between the x/y[/z] groups). The spacing may be non-uniform (the
grid is "stretched"), which the viewport's mesh overlay draws faithfully:

```
[mesh]
source = "file"
path = "grid.vtk"     # or "coords.txt"
```

The solver runs on this grid; the console's mesh tab renders imported
faces with the non-uniform spacing visible, over the field.

A clustered mesh can also be declared inline instead of imported: set
`mesh.source = "clustered"` with `first_cell` (the height of the first cell
at the refined side) and `growth` (the geometric ratio of each next cell),
optionally `cluster = "wall"`, `"channel"`, or `"top"` for which side is
refined. The flat-plate case under `cases/` uses this for its
boundary-layer grid.

Gmsh `.msh` (v2.2 ascii, triangles/quads) is also supported for 2D: the mesh
is parsed, built into a cell-centered unstructured grid, and solved with the
face-based `advance_ugrid` solver. `source = "file"` + `path = "mesh.msh"`
routes the case to that path automatically.

## Next

Read the [[architecture|architecture]] notes to understand the layers, then
follow the [[numerics/1d/euler|1D]] → [[numerics/2d/euler|2D]] →
[[numerics/3d/solver|3D]] solver pages in order. The
[[formats/case-toml|case format reference]] is the authoritative spec for the
case file.
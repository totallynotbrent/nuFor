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

| What you want | Key fields to change |
|---|---|
| A 2D/3D blast | `equations = "euler_2d"` or `"euler_3d"`, add `mesh.ny`/`nz`, `initial_condition.type = "blast"` |
| A custom mesh | `mesh.source = "file"` + `mesh.path` to a list of cell centers |
| A different end time | `time.final_time` |
| More resolution | `mesh.nx` (and `ny`/`nz`) |

The full schema is documented in [the case format reference](formats/case-toml.md).

## See it in the browser

```bash
cargo run --release --bin nufor -- serve 8060
```

Then open `http://localhost:8060`. The page is a single CFD workspace: a flow
viewport (2D animation or 3D slices with a toggleable mesh overlay), a left
sidebar for simulation and display controls, and a diagnostics dock. The **Run**
button solves your configured case and saves `results/<name>.vtk`.

## Bring your own mesh

For 2D/3D, set `mesh.source = "file"` with `mesh.path` pointing at a
**structured** mesh file — an ascii VTK `RECTILINEAR_GRID`/`STRUCTURED_POINTS`
block, or a plain `coordinates` file (one face coordinate per line, a `#` or
blank line between the x/y[/z] groups). The spacing may be non-uniform (the
grid is "stretched"), which the viewport's mesh overlay draws faithfully:

```
[mesh]
source = "file"
path = "grid.vtk"     # or "coords.txt"
```

The solver runs on this grid, and the web UI's "mesh file" field (Display →
mesh file → load mesh) renders the imported faces — non-uniform spacing
visible — over the field.

Gmsh `.msh` (v2.2 ascii, triangles/quads) is also supported for 2D: the mesh
is parsed, built into a cell-centered unstructured grid, and solved with the
face-based `advance_ugrid` solver. `source = "file"` + `path = "mesh.msh"`
routes the case to that path automatically.

## Next

Read the [architecture](architecture.md) notes to understand the layers, then
follow the [1D](numerics/1d/1d-euler.md) → [2D](numerics/2d/2d-euler.md) →
[3D](numerics/3d/3d-solver.md) solver pages in order. The
[case format reference](formats/case-toml.md) is the authoritative spec for the
case file.
---
title: Mesh diagnostics
---

nuFor validates a mesh before it is solved on. The diagnostics pass runs on the
unstructured `Ugrid` and reports the geometry the solver depends on: cell
volumes, how much those volumes vary, whether every cell closes under the
signed-face sum, whether any interior face points the wrong way, and whether
any cell is degenerate or negative.

## What is checked

`Ugrid::diagnostics()` returns a `MeshDiagnostics` summary with four groups:

- **Counts.** Number of cells, faces, interior (shared) faces, and boundary
  faces.
- **Area stats.** `min` / `mean` / `max` cell area and the `stretch` ratio
  (`max/min`). A stretch of 1.0 means every cell is the same size; a large
  stretch flags badly scaled cells.
- **Closure residual.** The signed face vectors of one cell should sum to zero.
  The pass reports the largest such residual over all cells (normalized by the
  cell perimeter) and how many cells are "open". A mesh with open cells is
  non-manifold or has inconsistently signed faces.
- **Face orientation.** An interior face points from its left cell to its right
  cell, so its area vector should agree with the center-to-center step. Faces
  that point the other way are flipped.
- **Negative volume.** Cells whose signed area is non-positive are counted.
  `from_cells` rejects these at build time; the diagnostics reports them anyway
  so a checked mesh is never silently solved on.

The summary carries a single `valid` verdict: true only when there are no
negative cells, no open cells, and no flipped faces.

## Using it

The quality pass runs at import time. From the command line:

```
nufor mesh-check case.msh
```

prints the table and a verdict, and exits non-zero for a mesh that has
problems:

```
cells: 4
faces: 12 (4 interior, 8 boundary)
area min/mean/max: 0.2500 / 0.2500 / 0.2500
area stretch (max/min): 1.00
negative cells: 0
largest closure residual: 0.00e0
open cells: 0
flipped interior faces: 0
verdict: valid mesh
```

The same summary is available to tooling as an HTTP route:

```
GET /api/mesh-stats?path=case.msh
```

## Why it matters

The finite-volume flux uses the face area vector and the cell volume as its
metrics. A degenerate or mis-oriented input produces a silent wrong answer or a
crash, not an error message. Running the diagnostics before solving turns that
into an explicit rejection, and the stretch and closure numbers give a first
read on mesh quality before a run starts.

## See also

- [[numerics/unstructured/unstructured-mesh||Unstructured mesh research]]
- [[numerics/unstructured/unstructured-fv||Unstructured finite volume]]
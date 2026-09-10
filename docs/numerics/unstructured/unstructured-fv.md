---
title: Unstructured finite-volume prototype
---

The research note laid out the destination: cell-centered finite volume on
arbitrary polygons, reusing the HLLC flux with a face-normal rotation. This
milestone is the working seed of that path -- a minimal mesh and a solver that
takes the flux you already know and points it at any face normal.

## The grid

`Ugrid` stores cells and faces; each face carries a signed area vector (its
length times the outward unit normal from the left cell toward the right cell),
which is the whole metric for the finite-volume flux. A boundary face has no
right neighbour and its area vector points out of its single cell. A
`cartesian_quads` constructor builds a quad mesh stored as arbitrary cells, so
the unstructured solver can be verified against the structured one.

## The step

`advance_ugrid` is a cell-centered finite-volume step. For every face it
rotates the two adjacent primitive states onto the face normal, applies the
existing HLLC flux in that rotated frame, and rotates the result back to global
components; the face flux then conservatively moves mass, momentum, and energy
from the left cell to the right. The time step is sized by a robust per-cell
length scale -- the inscribed width `area / (half perimeter)` -- because that is
what keeps an anisotropic cell (a tall thin quad) from being over-stepped by the
square-root-of-area measure the structured grid could afford.

## Verification

The Sod tube, run on a 100-by-4 quad mesh stored as arbitrary cells, reproduces
the classic solution: density spans 1.0 past 0.125, stays positive, and the
star plateau reads 0.418 against the exact 0.426 (first-order reconstruction,
so a little extra smearing is expected). The unstructured path is a front-end
onto the same physics, not a separate solver.

## See also

- [[numerics/unstructured/unstructured-mesh||Unstructured meshes]]
- [[numerics/unstructured/spalart-allmaras||Spalart-Allmaras]]
- [[numerics/2d/2d-euler||2D Euler solver]]

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
from the left cell to the right. The time step is sized by a per-cell
length scale -- the inscribed width `area / (half perimeter)` -- because that is
what keeps an anisotropic cell (a tall thin quad) from being over-stepped by the
square-root-of-area measure the structured grid could afford.

## Verification

The Sod tube, run on a 100-by-4 quad mesh stored as arbitrary cells, reproduces
the classic solution: density spans 1.0 past 0.125, stays positive, and the
star plateau reads 0.418 against the exact 0.426 (first-order reconstruction,
so a little extra smearing is expected). The unstructured path is a front-end
onto the same physics, not a separate solver.

A case whose `mesh.path` ends in `.msh` runs on this solver. A uniform initial
condition fills every cell and matches the inflow boundary, so an external
flow over an imported body works the same way the structured cases do. Both
closures run here: `eos = "eqair"` threads the equilibrium-air tables into
every cell and boundary ghost (a freestream pressure is inverted to its table
energy once per step), while `eos = "perfect"` stays on the gamma law. The
species lanes still have no unstructured march and are rejected loudly.

The reconstruction knob carries over: `reconstruction = "muscl"` builds
green-gauss gradients per cell and extrapolates face states with a two-cell
minmod limiter; `first_order` keeps the historical piecewise-constant states.

The gate for the closure plus reconstruction is a hypersonic ramp channel on
a 2,880-triangle mesh: Mach 8 air over a 14 degree ramp forms the oblique
shock, holds it steady across the last quarter of the run, and lands within
the coarse-mesh band of the theta-beta-M estimate (first order peaks at
3,530 Pa against roughly 3,900, muscl at 3,638). The limiter is not
transverse-invariant on triangles, so a small density undershoot rides the
shock foot (about 13 percent below freestream); the flow stays positive and
the shock position is unaffected.

The heartbeat's peak-pressure read used to carry an extra density factor, which
made every `p_max` the console ever showed for the unstructured runs
meaningless. It now matches the EOS closure exactly. Peak pressures from the
structured runs' own vtk fields were always computed through the model closure
and are unaffected.

## See also

- [[numerics/unstructured/unstructured-mesh|Unstructured meshes]]
- [[numerics/unstructured/boundary-groups|Boundary groups]]
- [[numerics/2d/sa-research|Spalart-Allmaras]]
- [[numerics/2d/euler|2D Euler solver]]

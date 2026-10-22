---
title: Boundary groups
---

Boundary faces now know what they are. A gmsh mesh carries physical groups on
its line elements, and the parser maps every boundary segment to the face that
shares its two nodes. Three groups drive the solver by convention:

- group 1, the wall. A slip wall takes the exact flux `[0, p·n, 0]` straight
  from the owning cell's pressure: no ghost state, no Riemann solve. A mirrored
  ghost state pushed through HLLC reads a few percent high at the wall because
  the sharpened wave speeds break the mirror cancellation, so the exact flux is
  both cheaper and correct.
- group 2, supersonic inflow. The face sees the freestream as its outside
  state; `ug.freestream` holds the case's initial condition.
- group 3, outflow. The outside state copies the owning cell, which is the
  supersonic outflow behavior.

Untagged boundary faces stay transmissive, so meshes from before the groups
keep working. `nufor mesh-check` prints the census: how many boundary faces each
group claimed and how many went untagged. A mesh whose wall loop has a gap
shows up there before it shows up as a leak in the run.

## Positivity floors

A first-order cell-average update can push a shock-layer cell slightly negative
at extreme Mach numbers. The march clamps every cell to a small fraction of the
freestream density and internal energy -- floors tied to the case's own initial
condition, set once at the start, never re-derived from the marching state. A
floor that recomputes from the current field ratchets downward: the wake
expands toward vacuum, the floor follows it down, and after a few thousand
steps a handful of near-empty cells carry ordinary pressure at density 1e-12.
Their sound speeds run to tens of kilometers per second and the global acoustic
time step collapses with them. Fixed floors tied to the freestream do not
ratchet.

## What the gates check

- Quiescent gas in a closed all-wall box holds exactly. A moving gas cannot
  stay uniform in a closed box -- it piles against the walls -- so the
  exact-hold invariant is for still air, where every flux row is zero at a wall
  and the pressure push cancels around each closed cell.
- A uniform wind in an open duct (inflow left, outflow right, slip walls top
  and bottom) holds exactly. This exercises all three groups at once with a
  nonzero normal velocity.
- Freestream inflow against a different interior pushes momentum into the
  domain, not out of it.
- The Sod tube on quad cells still lands on the classic plateau.

## See also

- [[numerics/unstructured/unstructured-fv|Unstructured finite-volume prototype]]
- [[numerics/unstructured/unstructured-mesh|Unstructured meshes]]

---
title: Supersonic cylinder
---

The supersonic cylinder is the first blunt-body case: a Mach-2 stream flowing
over a cylinder on the structured cartesian grid, producing a detached bow
shock that stands off ahead of the stagnation point. It exercises the
immersed-body path and the supersonic boundary conditions against a classic
high-speed flow.

## The immersed body

A circular body is cut out of the Cartesian mesh. Cells whose centres fall
inside the cylinder are frozen to a quiescent reference state, and the ring of
fluid cells adjacent to the surface has the velocity component normal to the
surface removed each step — a slip wall (no flow through the boundary). The
tangential velocity, density, and pressure are carried over, matching the
zero-gradient expectation at a solid wall. This lets the plain Cartesian
solver represent a curved wall without a body-fitted mesh.

## Working the flow

The domain runs a Mach-2 freestream in `+x` over a cylinder centred at
`(1.2, 0.5)` with radius `0.25`. Supersonic inflow on the west face, supersonic
outflow on the east, and slip walls on the top and bottom close the domain.
The solver runs to a quasi-steady state and the immersed body is imposed after
every step.

## What is verified

A detached bow shock forms ahead of the nose, and it is checked on the
stagnation line:

- **Upstream stays freestream.** Density ahead of the shock holds at `1.00`
  (the inflow is undisturbed).
- **A strong compression appears.** The stagnation-line density rises toward
  the normal-shock Rankine-Hugoniot ratio `rho2/rho1 ≈ 2.67` at Mach 2, with
  the stagnation-point peak a little above it (about `3.0`), matching the
  physical stagnation density behind the shock.
- **The shock is detached.** The compression front stands several cells ahead
  of the surface rather than attaching to the body, which is the signature of
  a blunt-body bow shock at supersonic speed.

The image below is the density field at a quasi-steady state: the dark disc is
the cylinder, the bright arc ahead of it is the bow shock, and the darker
region behind is the low-density wake.

## See also

- [[numerics/2d-boundaries|2D boundaries]]
- [[numerics/oblique-shock|Oblique shock]]
- [[numerics/shock-reflection|Shock reflection]]
- [[numerics/2d-euler|2D Euler solver]]
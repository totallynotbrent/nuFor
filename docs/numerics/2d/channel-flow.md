---
title: Laminar channel flow
---

With viscous terms in place the solver finally models a real viscous boundary
layer. The simplest such flow, steady incompressible pressure-driven flow
between two parallel no-slip plates, has an exact parabolic solution (Poiseuille
flow), so it is the natural first validation.

## The no-slip wall

A viscous flow against a solid wall must satisfy no slip: the fluid at the
wall moves with the wall, so its velocity vanishes there. The solver adds a
`NoSlipWall` boundary condition that mirrors both velocity components (so the
wall-average velocity is zero) while copying density and pressure. The shear
stress that develops near the wall is what Poiseuille's parabola comes from.

## Poiseuille's parabola

For a pressure gradient `dp/dx` (or its equivalent body force `B`), a no-slip
channel of height H carries a steady velocity

    u(y) = (B / 2 mu) y (H - y)

which peaks at `u_max = B H^2 / (8 mu)` mid-channel. The driving body force is
picked so the peak is small (well into the low-Mach, essentially incompressible
range), and the profile is entered as the initial condition.

## Verification

`tests/channel_flow.rs` lays the exact parabola into a channel clamped between
two no-slip walls, applies the matching body force, and marches the flow with
viscosity. Because the parabola is a steady solution of the discrete operator,
the solver must hold it in place -- and it does, to better than a percent over
the run, with density staying positive throughout. The profile error is
dominated by the accumulated drift of the body-force-driven, mildly
compressible flow, not by the viscous discretisation itself.

The same hold runs on a wall-clustered channel (`channel_grid2d`, first cell
hundreds of times finer than the mid-channel spacing): the position-aware
stencils keep the balance exact there too, which is what lets boundary-layer
cases put their resolution where the layer is.

## See also

- [[numerics/2d/viscous||Viscous terms]]
- [[numerics/2d/2d-euler||2D Euler solver]]
- [[numerics/foundations/positivity||Positivity]]

---
title: Blasius flat plate
---

The laminar flat plate is the first flow with an exact answer that exercises
everything the wall treatment promises: a no-slip wall, a clustered
wall-normal grid, a viscous layer that must neither thicken nor thin, and a
skin-friction coefficient with a closed form. The Blasius similarity
solution is that answer.

## The similarity solution

For a stream at speed `u_inf` over a plate starting at the leading edge, the
layer is self-similar in the variable

    eta = y sqrt(u_inf / (nu x))

with f(eta) obeying

    f''' + 0.5 f f'' = 0,  f(0) = f'(0) = 0,  f'(inf) = 1

The solver integrates this once at construction (fixed-step RK4 over
[0, 12] with the shooting value f''(0) = 0.3320573), stores the table, and
interpolates it monotonically. The table reproduces the textbook f' rows to
1e-4 and closes on f'(12) = 1 to roundoff, so the profile the boundary
condition prescribes is the published one, not a paraphrase of it.

The quantities that matter downstream:

- u = u_inf f'(eta)
- v = 0.5 u_inf (eta f' - f) / sqrt(Re_x)
- cf = 2 f''(0) / sqrt(Re_x) = 0.664 / sqrt(Re_x)
- delta99 = 5.0 x / sqrt(Re_x)

## The hold

`tests/blasius_hold.rs` lays the exact field into a wall-clustered grid
(leading edge upstream of the domain, inflow plane at Re_x = 40, about 27
cells across the layer at the exit), feeds the same analytic profile through
a ProfileInflow at the west face, and marches several flow-throughs.

The layer is held to about two percent of the freestream everywhere, and
the discrete wall shear matches the cf correlation to within a few percent
over the interior stations. The first stations sit in the inflow corner,
where the incoming profile and the wall meet, and the last ones sit in the
outflow corner; both carry the compatibility transient every bounded-domain
validation has, so the assertion covers the interior band where the physics
is clean.

## What it validates and what it does not

Validated: the wall treatment produces the right shear on a clustered grid,
the profile inflow feeds a layer consistently (the inviscid ghosts, the
reconstruction, and the viscous ghost gradients all use the same profile),
and the reported skin friction is the derivative the solver actually applied.

Not validated: turbulence. This is the laminar baseline. The SA model's
quantitative check at turbulent Reynolds numbers against the correlations
and the NASA TMR data is the remaining milestone, and this page is the
template for how that one will be judged.

## See also

- [[numerics/2d/2d-boundaries||2D boundary conditions]]
- [[numerics/2d/viscous||Viscous terms]]
- [[numerics/2d/channel-flow||Channel flow]]
- [[numerics/unstructured/sa-implementation||SA implementation]]

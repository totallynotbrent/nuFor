---
title: Turbulent flat plate
---

The turbulent flat plate closes the RANS story the Blasius page set up: the
Spalart-Allmaras model carrying a genuinely turbulent layer at Reynolds
numbers where the flat-plate correlations apply, judged against the same
published references production codes use.

## The operating point

The case runs at conditions chosen for the physics, not for prestige:
`Re_x` from about 5e4 to 5e5 along a 1.2 long plate, Mach 0.17, the SA
farfield value `nu_tilde_inf = 3 nu`, and a wall-clustered grid with a
first cell at `y+` of roughly 1 at the exit station. The virtual leading
edge sits 0.2 upstream of the inflow plane, because the layer is singular
at the edge itself; the inflow prescribes the analytic Blasius profile at
that station, so the plate enters the domain laminar and transitions
naturally on the wall.

## The march

The steady layer is reached by the time-accurate march, not by any
acceleration trick: about four flow-throughs (`t = 4`, on the order of an
hour at this resolution on a single core) brings the reported `cf(x)` to a
plateau, with `nu_tilde` growing to some hundreds of `nu` at mid-chord,
which is the signature of a resolved turbulent layer rather than a
laminarization.

An experimental local-time-stepping flag exists
(`numerics.local_time_stepping`): each cell advances at its own stability
bound instead of the global one. Its fixed-point algebra is verified (it
lands the same steady state as the march on the contract tests), but on
this stiff wall layer the iterator does not converge, so the march remains
the validated path and the flag stays documented as experimental.

## Judging the result

Two references bound the answer honestly:

- the Schlichting power law `cf = 0.0592 Re_x^(-1/5)` (valid roughly
  5e5 to 1e7, so it is stretched at the low end of this band), and
- the NASA TMR CFL3D SA solution, the same model solved by a production
  code on a far finer grid.

In the transition band the references disagree with each other by more than
the tolerance anyone would call agreement: CFL3D's own SA `cf` runs 12 to 15
percent below the power law at `Re_x` near 1e6, because the correlations
assume a turbulent layer from the edge while the SA solution earns its
transition on the wall. So the honest comparison is against the model
reference, with the correlations as context.

The measured mid-domain stations sit 12 to 19 percent above the power law
over `Re_x` from 1.3e5 to 3.3e5, which is the same sign and a similar
magnitude as the CFL3D-versus-correlation spread. The first and last few
stations sit in the inflow and outflow corners, where the bounded domain
forces compatibility the unbounded problem does not have; the assertion
band excludes them, as in the laminar hold.

## What it validates and what it does not

Validated: the SA transport sustains a turbulent layer with the published
constants and closures, the wall treatment holds up at `y+` near 1, and
the reported `cf` is the shear the solver actually applied at every
station.

Not validated: grid convergence of the turbulent solution (one grid, not
a family), the log-law `u+` profile in inner scaling (the TMR data for it
is banked, the comparison is future work), and Reynolds numbers above
about 5e5 where this grid loses resolution and the march becomes too
expensive to iterate on. Local time stepping does not yet converge on
stiff wall layers and is not on the validated path.

## See also

- [[numerics/2d/blasius||Blasius flat plate]]
- [[numerics/unstructured/sa-implementation||SA implementation]]
- [[numerics/2d/2d-boundaries||2D boundary conditions]]

---
title: Keeping the solution positive
---

A compressible solver that allows density or pressure to go negative is not a
solver — it is a crash with a stack trace. The 2D Euler module lives or dies by
whether the reconstruction, the wave speeds, and the time step can be trusted
to keep rho and p strictly above zero, so the machinery for that deserves its
own note.

## Why anything goes negative at all

The exact intercell average temperature of a compression is fine, but shocks
and near-vacuum states are where discrete schemes slip. Write a cell in
primitive variables and the flux is a smooth function of them; the problem is
that the reconstruction extrapolates those primitives -- and extrapolation, by
construction, is not bounded on the cell. A limiter that mis-ranks waves or a
CFL that is too optimistic can overshoot exactly the cell that was already the
densest or most rarefied, and the full multi-stage update then interpolates
through a non-physical, negative intermediate.

## The two floors that matter first

The scheme clamps rho and p to a small positive floor (about 1e-12) before they
touch a square root or a division, so a bad guess degrades to a slightly wrong
but still finite state rather than a NaN that poisons the whole field. That is
a safety net, not a fix -- the second, structural protection is the time step.

## CFL and positivity

The global time step is the product of a CFL number and the fastest wave speed
in the domain. For a centered scheme with a limiter this is a bona-fide
positivity constraint: too high a CFL ratio and the local Riemann solution
coalesces into a state that is no longer a convex combination of the two cell
averages, and positive preservation is lost. The code keeps the default CFL at
a conservative 0.5 for the multi-stage update, which is well inside the
a-priori bound for positivity of the HLLC-family fluxes, and the diagnostics
layer (`check_physical`) watches for the failure mode anyway, reporting a blow
up reason instead of propagating garbage.

## Where the linearization can mislead

Walking the shock relations in the verification suites, the single most common
surprise was a pressure estimated as slightly negative by an early Davis-style
wave-speed estimate while the true Riemann problem was perfectly fine. Picking
the pressure-based (PVRS) wave speeds -- a closed-form estimate of p* from the
left and right states rather than a crude bounding of the contact -- removed
the systematically wrong ranking that manufactured those negatives in the first
place. The floor stays for the corner cases; the wave-speed choice means it is
almost never reached.

The practical rule the code follows: preserve positivity structurally where you
can (wave speeds, CFL, limiting), and keep a floor + diagnostics for the
corners you cannot see coming.

## See also

- [[numerics/eos|Ideal-gas EOS]]
- [[numerics/flux-hll|HLL flux]]
- [[numerics/time-step|CFL time-step]]
- [[numerics/diagnostics|Diagnostics]]

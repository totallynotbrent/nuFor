---
title: HLL versus HLLC
---

The solver carries two flavours of Riemann flux: the classic HLL scheme used by
the 1D solver, and the HLLC scheme used by the 2D solver. They share the same
family tree and differ in one structural choice that decides how crisply they
resolve a contact surface.

## The two-wave and three-wave splits

HLL approximates the Riemann problem by two waves, a left wave with speed S_L
and a right wave with speed S_R, and collapses everything between them into one
constant intermediate state. That is the whole approximation: the structure
between the two waves -- in particular the contact discontinuity -- is smeared
away, its wave speed inferred from nothing but the two bounding states.

HLLC keeps the contact. It adds a third, middle wave (the contact, with speed
S*) and resolves genuinely distinct constant states on either side of it. The
cost is one extra star-state estimate per interface and a small amount of book
keeping; the benefit is that the contact and the two maze-like states behind it
are reproduced instead of levelled.

## What each buys

For shocks and supersonic regions the difference is small: both HLL and HLLC
capture strong discontinuities cleanly and stay positive for the same CFL
pool. The gap shows up at a contact -- the middle discontinuity in the Sod
tube, a material interface, a density front. There HLL rounds the plateau
(diffusive), while HLLC holds the jump to the grid resolution. On the Sod tube
the 2D HLLC path matches the reference solution to well under a percent where
the 1D HLL path shows a visibly broadened density plateau.

|            | HLL        | HLLC        |
|------------|-----------|-------------|
| waves kept | 2 (outer) | 3 (keeps contact) |
| contact resolution | smeared | crisp |
| cost per interface | cheapest | slightly more |
| positivity under CFL | yes | yes |
| the solver's use | 1D Euler | 2D Euler |

## The practical rule

Default to HLLC whenever contact/interface resolution matters -- which is most
real aerodynamics, and everything involving boundaries between species or
states. Reach for HLL when the per-interface cost is everything, the flow is
highly supersonic with no resolved contacts, or you want the cheapest possible
scheme for a code with a tight budget. The two live side by side in the solver
precisely so that cost and accuracy can be traded per use rather than fixed.

## See also

- [[numerics/flux-hll|HLL flux]]
- [[numerics/2d-euler|2D Euler solver]]
- [[numerics/1d-euler|1D Euler]]

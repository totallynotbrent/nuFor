---
title: 2D Euler solver
---

The two-dimensional extension of the finite-volume solver: an HLLC Riemann
flux with MUSCL reconstruction and a van Leer limiter, time-marched to second
order with an explicit two-stage (Heun) step. It builds directly on the 2D
state, grid, and EOS modules.

## Flux

`hllc_flux` computes the three-wave HLLC flux across an axis-aligned face from
the left/right primitive states. Wave speeds use a PVRS pressure estimate so
a rarefaction head is not under-ranked by a low-density neighbour, the star
state follows Toro's closed form, and the pressure/density are floored to a
small positive value so a reconstruction overshoot near a shock cannot poison
the sound speed. The flux is returned in global `(mass, mx, my, energy)`
components.

## Reconstruction

`face_states` reconstructs limited left/right interface values for each
primitive with the van Leer (harmonic-mean) limiter, which is TVD and keeps the
reconstructed value within the local data. With the limiter disabled the scheme
falls back to piecewise-constant, first-order HLLC — the two are separated in
the tests to show the order the limiter buys.

## Time stepping

`advance2d` is a single explicit step (CFL-limited on the fast-characteristic
speed over the cells); `advance2d_rk2` is the two-stage Heun predictor-corrector
that gives second order in time so the MUSCL spatial accuracy is not masked by a
first-order integrator.

## Verification

- Smooth isentropic vortex (Yee et al.) advected supersonically: the MUSCL
  scheme converges at about second order (`order` measured just above 2) while
  plain HLLC stays first order.
- A two-dimensional Sod tube aligned with x reproduces the one-dimensional
  solver to about one part in a hundred on every row.
- Mass and energy are conserved to machine precision and density stays
  positive.

The boundary handling here is simple transmissive ghost cells; supersonic
inflow/outflow, slip-wall, and symmetry conditions are a separate milestone.

## See also

- [[numerics/2d|2D foundations]]
- [[numerics/oblique-shock|Oblique shock]]
- [[numerics/hll-vs-hllc|HLL vs HLLC]]
- [[numerics/2d-boundaries|2D boundaries]]

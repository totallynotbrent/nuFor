# Sod shock tube

The classic one-dimensional Riemann problem: a high-pressure gas on the left
expanding into a low-density, low-pressure gas on the right, producing a left
rarefaction, a right-moving contact, and a right-moving shock.

## Setup

- Domain [-0.5, 0.5], gamma = 1.4.
- Left: rho = 1.0, u = 0, p = 1.0.  Right: rho = 0.125, u = 0, p = 0.1.
- Time: t = 0.22, transmissive outer boundaries.

## Verification

The scheme is compared against the exact solution of the 1D Riemann problem
(see the `exact` module). The exact star pressure and contact velocity match
the published values (p* = 0.30313, u* = 0.92745), and the numerical density
field agrees with the exact solution with an L1 error in the low one-percent
range on a 400-cell mesh. The L1 density error shrinks as the mesh refines,
confirming first-order convergence.

The Sod test rides on the solver, flux, and CFL tests already landed, and is
the first rung of the verification ladder.

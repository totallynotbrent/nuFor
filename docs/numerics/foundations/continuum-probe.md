---
title: Continuum-validity probe
---

Rung 7: the run-log line that answers whether the Navier-Stokes
closure the solver integrates is even valid where it is being
integrated. I added it because every rung on the ladder marches a
continuum solution, and the ladder's later rungs (catalytic walls,
combustion) inherit that assumption; the probe states it explicitly
per run instead of leaving it implicit.

## Bird's breakdown parameter

The continuum closure assumes local near-equilibrium. Bird's
criterion makes that a number: instead of the global
Knudsen number on a body scale, it uses the local
gradient length as the macroscopic scale, so it can flag a
shock foot or a sharp leading edge even when the global
Kn says continuum:

$$\mathrm{Kn}_{GLL}(q) = \frac{\lambda}{q}\,|\nabla q|$$

with $q$ the density (Bird's original form) or the temperature,
and $\lambda$ the local mean free path. Boyd et al. (Physics of
Fluids 7, 1995) and Lofthouse and Boyd (AIAA 2007-3905) use
$\mathrm{Kn}_{GLL} > 0.05$ as the continuum breakdown threshold
when they switch between CFD and DSMC domains in hybrid solvers;
I report the same threshold so a nuFor run states where its
continuum assumption holds and where it does not.

## Mean free path

The probe needs $\lambda$ per cell, which means it needs collision
data or a shortcut. I use the shortcut: for a hard-sphere gas the
viscosity is $\mu = \rho\,\lambda\,\bar{c}/2$ with the mean
molecular speed $\bar{c} = \sqrt{8RT/\pi}$, so

$$\lambda = \frac{2\mu}{\rho\,\sqrt{8RT/\pi}}$$

which stays well-defined for any gas whose viscosity law is known
(Sutherland here) and needs no collision diameters. The velocity is
in the viscosity relation, so $\lambda$ carries the same
temperature dependence the frozen-Pr analysis on the Fay-Riddell
page used, and the two diagnostics agree with each other.

## What the run log prints

After every 2D or axisymmetric run the log carries:

```
continuum probe (bird gradient-length knudsen):
  kn_inf (lambda_inf/l_ref) 1.215e-7
  max kn_gll(rho) 3.206e-2 at (1.925, 0.825)
  max kn_gll(t) 2.623e-3
  cells above breakdown threshold 0.0 pct
```

- `kn_inf` is the freestream mean free path over the case's
  reference length, the global number.
- the two `max kn_gll` lines are the worst cell of each form,
  with its location.
- the percentage is the fraction of flow cells above 0.05 in
  either form, the field-wide statement of the same criterion.

Solid cells are skipped, and the gradient stencil steps over them:
a masked body's pinned interior values never reach the quotient.
The reference length comes from the case's `[physics.reference] l`
when it is set, else the domain's x extent.

## Reading it

The Casseau and capsule runs sit far below the threshold
throughout ($\mathrm{Kn}_{GLL} \sim 10^{-2}$ at the worst shock
cells, $\mathrm{Kn}_\infty$ six to seven orders lower), which is
the expected verdict for those densities and no surprise. The
probe earns its keep at the top of the atmosphere: the trajectory
sweep's highest-altitude station is where the stagnation region
first approaches the 0.05 line, and a reentry case run at skip-
entry densities will state exactly which cells crossed it. That is
the input the DSMC rung would need (rung 7 exists to tell the
hybrid story honestly, and the threshold is the industry-standard
0.05 so the numbers here compare directly with the hybrid-CFD
literature).

## Module

`nufor-core::continuum_probe`. Unit tests cover the uniform-field
zero, a linear ramp against the hand value, solid-neighbor
skipping, and the mean-free-path scaling laws.

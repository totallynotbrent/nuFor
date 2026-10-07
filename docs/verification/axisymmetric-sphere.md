---
title: Sphere validation
---

The axisymmetric solver's contract case: Mach-2 flow over a sphere,
checked against the Ambrosio-Wortman bow-shock standoff correlation and
the Rankine-Hugoniot post-shock states. It is the case that separates
the annular solver from the planar one, because the standoff of a sphere
is about twice that of a cylinder and only the geometric source terms
make it so.

## Setup

A sphere of radius 1 centred at (2, 0), masked onto the grid with the
rung-2 solid-cell mask. The half-domain is $y \in [0, 3]$ with the
symmetry axis on the south face (`SlipWall`, which is the axis
condition: the ghost's negated radial velocity makes the axis flux
mass-free). Mach-2 inflow on the west face, outflow on the east. The
grid is 180x90, and the march runs to $t = 4$.

## What is checked

- Pre-shock density stays at 1.000: the inflow reaches the shock
  undisturbed, so nothing upstream of the bow shock has leaked.
- The post-shock density peaks at 3.03, inside the Rankine-Hugoniot
  band: 2.667 is the normal-shock jump at Mach 2, and the stagnation
  rise puts the peak a little above it.
- The standoff measures 0.317, within the 40 percent tolerance of the
  Ambrosio-Wortman reference 0.321. The coarse grid and the staircase
  sphere shift the effective nose, which is why the tolerance is loose:
  at this grid the measured miss is 1.2 percent.

## The reference

The Ambrosio-Wortman correlation is a real 1962 result (ARS J. 32(2)
281, "Stagnation-Point Shock-Detachment Distance for Flow Around
Spheres and Cylinders in Air"), with the sphere form

$$
\delta/R = 0.143 \exp(3.24/M^2).
$$

It was pulled from the literature during the sprint, not recalled, and
the numbers above are checked against its Mach-2 value rather than
against a remembered constant.

## See also

- [Axisymmetric Euler](../numerics/2d/axisymmetric.md)
- [Regression suite](regression-suite.md)

---
title: 3D channel validation
---

The 3D Poiseuille hold: the analytic parabolic profile between two
no-slip plates, carried in a three-dimensional slab, driven by the
equivalent pressure-gradient body force, and verified steady.

## Setup

A slab of channel $[0, 0.5] \times [0, 1] \times [0, 0.25]$: no-slip
plates at $y = 0$ and $y = 1$, transmissive everywhere else. The
initial state is the exact parabola $u(y) = 4 u_{max} y (1 - y)$ with
$\mu = 0.05$, $u_{max} = 0.05$, and the body force $B = 8 \mu u_{max}$
that balances $\mu\, \partial^2 u / \partial y^2$ everywhere. The march
runs `advance3d_visc_rk2` for a physical time of 0.05, about ten steps
at the diffusion-capped increment.

## What is checked

The recovered $u(y)$ at the streamwise mid-plane must remain the
parabola at every spanwise station. The profile must neither decay
(the viscous drain unbalanced) nor grow (the body force unbalanced),
and must not vary in $z$ (the operator must not manufacture spanwise
structure from a spanwise-invariant field). Density must stay positive.

Measured: the worst profile error over the whole cross-section is about
$10^{-3}$, two percent of $u_{max}$, and spanwise invariance holds to
roundoff.

## See also

- [3D Navier-Stokes](../numerics/3d/3d-viscous.md)
- [Channel flow (2D)](../numerics/2d/channel-flow.md), the 2D original

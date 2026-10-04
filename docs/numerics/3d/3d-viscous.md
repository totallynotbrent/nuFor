---
title: 3D Navier-Stokes
---

The 2D viscous operator grown into three directions: the full
Navier-Stokes diffusive flux on the uniform structured 3D grid, with
no-slip walls on any of the six domain faces.

## The operator

`add_viscous3d` computes the diffusive fluxes of the three momentum
components and energy across every interior face. Each face family
carries its stress row with the Stokes $-2/3 \mu \nabla \cdot u$ term:
for an x-face, $\tau_{xx} = \mu(2 u_x - \tfrac{2}{3}\nabla\cdot u)$,
$\tau_{xy} = \mu(u_y + v_x)$, $\tau_{xz} = \mu(u_z + w_x)$, and the
energy flux is the face velocity dotted into that row plus heat
conduction $\kappa \partial_x T$. Gradients are cell-centered and
averaged to faces; the Prandtl number links conduction to viscosity.

## Walls on all six faces

`Walls3d` marks any subset of the six domain faces solid; every open
face stays transmissive, matching the inviscid step's ghost treatment.
Wall faces carry the one-sided quadratic-consistent shear
$\partial p/\partial n = (9 p_0 - p_1)/(3 \Delta)$, exact for profiles
quadratic near the wall, so the wall flux and the wall-cell balance
agree the same way the 2D operator guarantees. Walls are adiabatic: no
heat flux and no wall-parallel velocity, so the energy flux through a
wall face vanishes. The sign convention follows the divergence form:
the low-side wall face enters its cell as a loss, the high-side wall
face enters with the physically negated tangential stress.

## The march

`advance3d_visc_rk2` is the 3D counterpart of `advance2d_visc_rk2`:
Heun's two-stage step, each stage an inviscid sweep capped by the
explicit diffusion bound followed by the viscous add, averaged at the
end.

## Validation

Two flows close the operator against analytic solutions, both using the
hold style: initialize the exact steady state and verify the march
keeps it rather than developing toward it.

- [Channel flow (3D)](../verification/viscous3d-channel.md): the
  Poiseuille parabola between two plates, held across the spanwise
  direction.
- [Square duct](../verification/viscous3d-duct.md): the
  double-Fourier series solution, held on all four lateral walls.

A third contract guards the operator's axis symmetry: a shear field
between y-walls drains x-momentum identically to the same field between
z-walls, which the channel case (a z-invariant flow) cannot see.

## See also

- [3D solver](3d-solver.md), the inviscid step this operator augments
- [Viscous terms](../2d/viscous.md), the 2D original

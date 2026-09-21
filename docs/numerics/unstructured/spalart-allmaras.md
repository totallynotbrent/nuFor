---
title: Spalart-Allmaras turbulence research note
---

## Question

The 2D and 3D Euler solvers are now viscous (the Navier-Stokes path). At
engineering Reynolds numbers a direct simulation is unaffordable, so the next
upgrade is a turbulence model. This note asks whether Spalart-Allmaras (SA) is
the right first RANS model for nuFor, and records what implementing it honestly
would require.

## Scope

A research note, not an implementation: the closure's assumptions, its single
extra field, the discretisation it needs, and the validation it should be
measured against. The decision is a roadmap statement, not a code commit.

## Findings

Spalart-Allmaras is a one-equation eddy-viscosity model: one PDE for a modified
turbulent viscosity that a transport equation pushes toward an algebraic eddy
viscosity far from walls. Its appeal for an open-source CPU solver:

- **one extra scalar** -- a single transported field (the modified viscosity)
  on top of the five conserved variables, so it is cheap and the existing
  finite-volume machinery carries it,
- **simple wall treatment** -- the no-slip walls already in the viscous solver,
  with an SA wall boundary that zeroes the field there,
- **a well-trodden reference implementation** -- the 1994 NASA-Ames baseline has
  a canonical form and published calibration constants, so a from-scratch
  version has an external answer to check against,
- **designed for aeronautics** -- it was invented for the transonic wing and
  airfoil problems nuFor's mesh path is pointed at.

The model's single PDE transports a modified turbulence variable with
production, destruction, and two non-linear wall-damping terms that are the
hard parts to get stable on a coarse mesh. The alternative two-equation models
(k-epsilon, k-omega, SST) add a second equation and model constants that are
arguably harder to calibrate correctly than SA's single field.

## Implementation sketch

An SA build slots into the existing viscous 2D/3D solver as:

1. a sixth conserved field `nu_tilde` carried by the same HLLC-style advection
   and the same viscous machinery,
2. a source term (production/destruction/wall damping) evaluated at each cell
   from the velocity gradients the solver already computes,
3. the eddy viscosity `mu_t = rho * nu_tilde * f_v1` added to the molecular
   viscosity in the momentum and energy viscous fluxes,
4. a no-slip wall boundary that sets `nu_tilde = 0`.

Case 9 (flat-plate boundary layer) from the SA-typical validation set is the
entry test: skin friction along a flat plate must track the published
correlation.

## Decision

Adopt Spalart-Allmaras as the first turbulence model when the turbulence
foundation milestone is built, ahead of two-equation models. It is the cheapest,
best-documented one-equation model and its reference implementation makes honest
verification tractable. It does not retroactively change the current
Navier-Stokes solver; it is the next additive layer.

## Verification implications

- SA must reduce to laminar flow when the transported viscosity is zero.
- The flat-plate boundary layer must reproduce the published skin-friction
  correlation, not just a plausible shape.
- The model must not turn a clean laminar channel (the existing Poiseuille
  test) into a turbulent one at the low Reynolds number where it is genuinely
  laminar -- a guard against the model being numerically over-active.

## Open questions

- Wall-damping near the no-slip face needs the distance to the wall; for the
  cell-centered Cartesian path that is a geometric lookup the mesh provides.
- Explicit time stepping with SA production can be stiff; the diffusive-time
  bound that already governs the viscous term likely governs this one too.

## Date

2026-06-10

## See also

- [[numerics/unstructured/sa-implementation||SA implementation]]
- [[numerics/unstructured/unstructured-fv||Unstructured FV]]
- [[numerics/2d/viscous||Viscous terms]]
- [[numerics/unstructured/unstructured-mesh||Unstructured meshes]]

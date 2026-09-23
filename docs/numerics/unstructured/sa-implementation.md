---
title: SA implementation
---

The Spalart-Allmaras model as it is implemented in nuFor: what runs per step,
how it couples to the mean flow, and the guards that keep an explicit scheme
stable.

## What advances per step

The state carries one extra scalar per cell, the working variable
`nu_tilde`. A `rans_2d_sa` case advances it together with the mean flow in a
second-order Heun step:

1. the inviscid Euler update runs with the usual HLLC flux and MUSCL slopes,
2. the viscous fluxes are added with the molecular viscosity plus the eddy
   viscosity `mu_t = rho * nu_tilde * fv1(chi)`, with the turbulent Prandtl
   number carrying the eddy part of the heat flux,
3. `nu_tilde` itself is transported: upwind advection on the mean-flow face
   velocities, central diffusion with the `(nu + nu_tilde)/sigma` diffusivity
   plus the `c_b2/sigma |grad nu_tilde|^2` cross term, and the
   production/destruction source evaluated from the velocity gradients the
   viscous pass already has.

Stage two re-evaluates all of it at the stage-one state and the Heun average
closes the step, so the coupling is consistent to the same order as the
laminar stepper.

## Closures

The closures follow the 1994 baseline without the trip terms, with the
standard guards the reference implementations use:

- `chi = nu_tilde / nu`, `fv1 = chi^3 / (chi^3 + c_v1^3)`
- `fv2 = 1 - chi / (1 + chi fv1)`
- `stilde = s + fv2 nu_tilde / (kappa^2 d^2)`, floored at `cs * s` so the
  freestream region where `fv2 < 0` cannot drive it negative
- `r = nu_tilde / (stilde kappa^2 d^2)`, clipped at 10
- `g = r + c_w2 (r^6 - r)`, `fw = g ((1 + c_w3^6)/(g^6 + c_w3^6))^(1/6)`
- production `c_b1 stilde nu_tilde`, destruction `c_w1 fw nu_tilde^2 / d^2`

Constants: `c_b1 = 0.1355`, `c_b2 = 0.622`, `sigma = 2/3`, `kappa = 0.41`,
`c_w1 = c_b1/kappa^2 + (1 + c_b2)/sigma`, `c_w2 = 0.3`, `c_w3 = 2`,
`c_v1 = 7.1`, `cs = 0.3`. The vorticity `s` in 2D is `|dv/dx - du/dy|`, the
planar magnitude of `sqrt(2 omega_ij omega_ij)`.

## Boundaries and walls

A solid side of the domain is a no-slip wall for the mean flow and zeroes
`nu_tilde` in the ghost layer, which is the SA wall condition. Open sides
copy the interior. The wall distance `d` is the perpendicular distance to the
nearest solid plane of the structured grid; cells far from any wall get the
domain size, where destruction vanishes as `1/d^2`.

## Stability under explicit stepping

The transported equation is stiff in two places, both handled:

- the source is bounded by construction: production is linear in `nu_tilde`
  and destruction quadratic, so `nu_tilde = 0` is an exact fixed point and
  the update never invents turbulence from nothing,
- the time step carries the SA diffusion bound
  `dt < 0.25 dx^2 / ((nu + nu_tilde)/sigma)`, the same shape as the viscous
  bound, and `nu_tilde` is clipped non-negative after every step.

The freestream has a small slow decay (at `chi = 3`, `fv2 < 0` floors
`stilde` and destruction slightly wins); this is the model's known
freestream behavior, not an instability.

## Verification

- every closure is unit-tested against hand-checked anchors: `fv1(c_v1) =
  0.5`, `g(0) = 0`, `fw(0) = 0`, `fw(1) = (65/66)^(1/6)`, the composite
  `c_w1`
- `nu_tilde = 0` reproduces the laminar viscous stepper exactly, so existing
  laminar cases are provably untouched
- a uniform freestream is a fixed point up to the known slow decay
- a Gaussian peak decays at the analytic heat-kernel rate
- a no-slip channel runs hundreds of coupled steps with the state physical
  and `nu_tilde` bounded
- the no-slip wall treatment carries the one-sided quadratic-consistent wall
  shear, exact for the linear sublayer and the laminar parabola: the
  Poiseuille hold balances to machine precision on the discrete stencil
- a flat-plate run develops a positive skin-friction distribution along the
  wall, reported per station by `nufor run` (cf rows after the summary)

The flat-plate *quantitative* validation (Cf(x) against the flat-plate
correlations at turbulent Reynolds numbers) still needs an inflow that
carries a boundary-layer profile; the uniform inflow state does not support
that yet. Wall-normal clustering is in: the solvers are metric-aware on
non-uniform rectilinear grids (stretched or clustered meshes solve with the
exact position-aware stencils, and the flat-plate case now ships with a
clustered mesh). Until the profile inflow lands, the plate case stays
qualitative: the layer grows, the reported Cf sits in a plausible band, and
the SA field stays bounded.

## See also

- [[numerics/unstructured/spalart-allmaras||The SA research note]]
- [[numerics/2d/viscous||Viscous terms]]
- [[formats/case-toml||The case file]]
---
title: Spalart-Allmaras model
---

The Spalart-Allmaras model as it is implemented in nuFor: what runs per step,
how it couples to the mean flow, and the guards that keep an explicit scheme
stable.

## What advances per step

The state carries one extra scalar per cell, the working variable
$\tilde\nu$. A `rans_2d_sa` case advances it together with the mean flow in a
second-order Heun step. The transported equation is the SA transport equation,

$$\frac{D\tilde\nu}{Dt} = c_{b1}\,\tilde S\,\tilde\nu + \frac{1}{\sigma}\nabla\!\cdot\!\left[(\nu + \tilde\nu)\nabla\tilde\nu\right] + \frac{c_{b2}}{\sigma}\left|\nabla\tilde\nu\right|^2 - c_{w1} f_w \left(\frac{\tilde\nu}{d}\right)^{2}$$

production, diffusion, cross-diffusion, destruction, and the eddy viscosity
closes on the mean flow as

$$\mu_t = \rho\,\tilde\nu\, f_{v1}(\chi), \qquad \chi = \frac{\tilde\nu}{\nu}$$

1. the inviscid Euler update runs with the usual HLLC flux and MUSCL slopes,
2. the viscous fluxes are added with the molecular viscosity plus the eddy
   viscosity above, with the turbulent Prandtl number carrying the eddy part
   of the heat flux,
3. `nu_tilde` itself is transported: upwind advection on the mean-flow face
   velocities, central diffusion with the $(\nu + \tilde\nu)/\sigma$
   diffusivity plus the $c_{b2}/\sigma\,|\nabla\tilde\nu|^2$ cross term,
   and the production/destruction source evaluated from the velocity
   gradients the viscous pass already has.

Stage two re-evaluates all of it at the stage-one state and the Heun average
closes the step, so the coupling is consistent to the same order as the
laminar stepper.

## Closures

The closures follow the 1994 baseline without the trip terms, with the
standard guards the reference implementations use:

$$f_{v1} = \frac{\chi^3}{\chi^3 + c_{v1}^3}, \qquad f_{v2} = 1 - \frac{\chi}{1 + \chi f_{v1}}, \qquad \tilde S = s + \frac{f_{v2}\,\tilde\nu}{\kappa^2 d^2}$$

with $\tilde S$ floored at $c_s\,s$ so the freestream region where
$f_{v2} < 0$ cannot drive it negative, and

$$g = r + c_{w2}(r^6 - r), \quad f_w = g\left(\frac{1 + c_{w3}^6}{g^6 + c_{w3}^6}\right)^{1/6}, \quad r = \frac{\tilde\nu}{\tilde S\,\kappa^2 d^2}$$

with $r$ clipped at 10. Production is $c_{b1}\tilde S\tilde\nu$ and
destruction $c_{w1} f_w \tilde\nu^2/d^2$.

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
  0.5`, `g(0) = 0`, `fw(0) = 0`, `fw(1) = 1` exactly (the log-layer
  equilibrium of the published damping function), and the composite
  `c_w1` constant against its definition
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

The turbulent flat-plate validation is done, and it sharpened the model
along the way: the `fw` damping function carried a stray `+1` in its
denominator that a constant-by-constant audit against the published model
caught, the Blasius inflow now refuses to anchor at the leading edge
where the similarity solution is singular, and the shipped plate case
moves its virtual leading edge upstream accordingly. The turbulent layer
is sustained at Reynolds numbers where the flat-plate correlations apply:
`nu_tilde` grows to hundreds of `nu`, and the discrete Cf sits 12 to 19
percent above the power law over the interior stations, which is the same
sign and scale as the spread between the CFL3D SA reference and the
correlations themselves. See [[numerics/2d/turbulent-plate|the turbulent
plate page]] for the full judgment, including what remains open (grid
convergence, the inner-scaling log-law profile, Reynolds numbers past
about 5e5).

## See also

- [[numerics/2d/sa-research|The SA research note]]
- [[numerics/2d/viscous|Viscous terms]]
- [[formats/case-toml|The case file]]
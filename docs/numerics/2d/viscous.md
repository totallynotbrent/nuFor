---
title: Viscous terms
---

The Euler equations are the inviscid limit of the Navier-Stokes equations.
Adding viscosity and heat conduction turns the sharp, non-diffusive solver into
a full compressible Navier-Stokes one, at the cost of second-order spatial
derivatives and a more restrictive time step. The viscous terms are layered on
top of the existing 2D Euler update rather than mixed into it, so the inviscid
legacy stays intact and a single flag turns viscosity on.

## The viscous flux

For a Newtonian fluid with a constant dynamic viscosity mu and a Prandtl
number Pr, the diffusive part of the conservation law is

    d/dt(rho)                = 0
    d/dt(rho u)              = div(tau)
    d/dt(rho e)              = div(u . tau) - div(q)

where the shear stress and heat flux are

    tau_ij = mu ( du_i/dx_j + du_j/dx_i - (2/3) div(u) delta_ij )
    q      = -k grad(T),       k = mu * c_p / Pr

and c_p = gamma/(gamma-1) (with R = 1) sets the conductivity from the viscosity
and the Prandtl number. Temperature is recovered from the ideal-gas relation
T = p/rho.

## Numerical treatment

`viscous2d.rs` computes cell-centred velocity and temperature gradients by
central differencing (one ghost layer mirrored at the boundary), averages them
onto each cell face, and the flux divergence is added to the same conservative
update as the inviscid flux, but with the sign a diverging diffusive flux
requires. The domain-edge faces carry zero viscous flux, so the viscous terms
conserve mass and total energy exactly (no transport leaks through the sides).

On a non-uniform rectilinear grid the same machinery runs metric-aware: each
cell takes the exact three-point Lagrange derivative at its true center over
its neighbors (quadratic-exact on any spacing, the plain central difference
when the axis is uniform), face values interpolate to the true face
position, and the divergence divides by the cell's own width. The wall
treatment generalizes the same way: the one-sided quadratic through the wall
value and the first two cell centers is evaluated on the true positions, so
the stretched Poiseuille hold balances to machine precision. Uniform grids
keep the original scalar fast path bit for bit.

A profile inflow side feeds the layer through this machinery too: the ghost
velocities come from the profile at each row's true height (the Blasius
profile or a loaded table), so the centred wall-normal gradient at the
inflow column sees the shear the incoming layer actually carries. That is
what keeps the Blasius hold balanced from the first cell.

## Stability

Explicit diffusion is unstable unless the time step respects the diffusive
bound, dt <= h^2 / (2 D), where the binding diffusivity D = mu*kappa/rho with
kappa = gamma/((gamma-1) Pr). Because heat conduction usually diffuses several
times faster than momentum (kappa ~ 4.9 for air), that conduction limit is the
one that binds. The viscous time stepper (`advance2d_visc_rk2`) clamps its
step to that limit with a safety factor, so turning on viscosity never silently
outruns the scheme's stability.

## Sutherland's law

Constant viscosity is the simple choice for a first viscous milestone, but a
temperature-dependent law is provided too:

    mu(T) = mu0 (T/T0)^1.5 (T0 + S) / (T + S)

with the air constants mu0 = 1.716e-5, T0 = 273.15 K, S = 110.4 K. It lowers the
computational cost of the channel-flow milestone's high-temperature cases.

## See also

- [[numerics/2d/channel-flow||Channel flow]]
- [[numerics/2d/2d-euler||2D Euler solver]]
- [[numerics/foundations/positivity||Positivity]]

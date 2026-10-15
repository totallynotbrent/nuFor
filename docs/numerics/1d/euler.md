---
title: 1D Euler
---

The first serious numerical milestone: a conservative, finite-volume 1D Euler
solver for an ideal gas, evolving density, momentum, and total energy and
recovering physical shocks, contacts, and expansions.

The solver integrates the 1D compressible Euler equations,

$$\frac{\partial}{\partial t}\begin{bmatrix} \rho \\ \rho u \\ E \end{bmatrix} + \frac{\partial}{\partial x}\begin{bmatrix} \rho u \\ \rho u^2 + p \\ u(E + p) \end{bmatrix} = 0, \qquad p = (\gamma - 1)\left(E - \tfrac{1}{2}\rho u^2\right)$$

in conservative finite-volume form: each cell stores cell averages, and the
update is the divergence of numerical face fluxes,

$$U_i^{n+1} = U_i^n - \frac{\Delta t}{\Delta x}\left(F_{i+1/2} - F_{i-1/2}\right)$$

where each face flux $F$ comes from the HLL solver between the reconstructed
left and right states, and the step $\Delta t$ obeys the CFL condition with
Davis wave-speed estimates.

## What is in place

- Conservative and primitive state with the conversion between them.
- A uniform 1D control-volume grid.
- The ideal-gas equation of state with pressure recovery and validity checks.
- The HLL flux with Davis wave-speed estimates, and CFL time-step control.
- The conservative finite-volume solver itself: it assembles ghost-cell states
  at the two boundaries (transmissive or reflective), forms the HLL flux at
  every face, advances the state with the explicit CFL time step, and reports
  the largest absolute conserved-variable change as the residual. It can log
  every step (step, time, dt, residual) and run until the residual drops below
  a tolerance, a time limit elapses, or a step budget is used up.

The solver reuses the Fortran kernels for the CFL step, state conversion, and
HLL flux; the boundary logic, the conservative update, the residual, and the
logging live in the Rust boundary.

## Verification

Covered by the euler test suite: a constant state stays at rest, reflective
walls conserve total mass and energy to rounding, a Sod shock tube stays
density-positive and evolves physical waves, and the first step uses exactly
the CFL time step.

## Verification ladder

Constant state → uniform advection/contact → Sod tube → Lax tube → stationary
shock → isentropic expansion. Each feature ships with a focused verification
case and evidence; a successful run or a screenshot is not verification.

## See also

- [[numerics/foundations/state-grid|State and grid]]
- [[numerics/fluxes/hll|HLL flux]]
- [[numerics/1d/sod|Sod verification]]
- [[verification/index|Verification ladder]]

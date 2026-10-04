---
title: Square duct validation
---

Laminar pressure-driven flow in a square duct, the 3D viscous
operator's strongest analytic test, because all four lateral walls
interact and no direction is privileged.

## The analytic solution

The steady solution of $\mu \nabla^2 u = -G$ on the square
$[-a, a]^2$ with $u = 0$ on the walls is the double-Fourier series

$$
u(y, z) = \frac{16 G}{\mu \pi^2}
  \sum_{m,n \text{ odd}}
  \frac{(-1)^{(m+n)/2}}{mn (k_m^2 + k_n^2)}
  \cos(k_m y) \cos(k_n z),
\qquad k_m = \frac{m \pi}{2a}.
$$

Two details of this reference deserve their own warning, both learned
by getting them wrong first: the wavenumbers are $m\pi/2a$ (the walls
sit at $\pm a$, not $\pm a/2$), and the prefactor is
$16 G / (\mu \pi^2)$, not $16 a^2 G / (\mu \pi^3)$. With both wrong the
series is still a plausible-looking smooth bump that is NOT the duct
solution: the discrete Laplacian of the field disagrees with the
analytic $\nabla^2 u = -G/\mu$ by 100 percent, and a hold test against
it fails in a way that looks like an operator bug.

The peak velocity checks against the literature value
$u_{max} = 0.295\, G a^2 / \mu$ to under one percent.

## The hold test

The duct $[0, 0.5] \times [-0.5, 0.5]^2$ is initialized with the
series evaluated at every cell center (40x40 modes), no-slip on the
four lateral faces, driven by the body force $G$, and marched with
`advance3d_visc_rk2` for a short physical time.

Measured: the worst profile error over the cross-section is 0.8 percent
of $u_{max}$, and the mean-to-max velocity ratio matches the analytic
field's own cell-sampled ratio to 0.1 percent. The ratio comparison is
made against the series sampled at cell centers, not the textbook
integral value 0.5699: center sampling shifts the ratio, and comparing
the simulation against the integral would measure the sampling, not
the solver.

The duct also motivated a companion contract: the operator must be
symmetric under the y-z axis swap, checked by draining a linear shear
field between y-walls against the same field between z-walls. The
channel case is z-invariant and cannot see a broken z-face momentum
path; the duct and this symmetry contract can.

## See also

- [3D Navier-Stokes](../numerics/3d/3d-viscous.md)
- [3D channel validation](viscous3d-channel.md)

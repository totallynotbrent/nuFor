---
title: Axisymmetric Euler
---

The 2D march revolved about an axis: every body of revolution (sphere,
sphere-cone, capsule) becomes a cheap, high-resolution 2D problem. The
module is `crates/nufor-core/src/axi.rs`, with `advance2d_axi` the single
step and `advance2d_axi_rk2` the Heun two-stage wrapper, composed over it
the same way `advance2d_rk2` composes over `advance2d`.

## The equations

In cylindrical coordinates (x axial, r radial) the axisymmetric Euler
equations keep the planar face fluxes unchanged, per unit area. What
changes is the finite volume: a planar cell is a rectangle, an
axisymmetric cell is an annulus, so the control volume shrinks toward
the axis and the radial flux difference must carry the face radii. In
divergence form,

$$
\partial_t U = -\partial_x F_x - \frac{1}{r} \partial_r (r F_r) + S,
\qquad S = (0, 0, p/r, 0),
$$

with $S$ the geometric source. Only the radial momentum component feels
it: the outer cylindrical wall of the annulus has more area than the
inner one, so even a pressure field uniform in $r$ pushes the ring
outward. In the discrete update on a cell with center radius $r_c$ and
radial face radii $r_-$ and $r_+$, the radial flux difference becomes
$(r_+ F[j+1] - r_- F[j]) / (r_c \, dr)$ and the radial momentum gains
the source $p / r_c$. The march reuses the planar flux machinery
unchanged and rewrites only the update.

## The axis boundary

The south face of the domain is the symmetry axis, and `Bc2d::SlipWall`
there is exactly the axis condition. The slip-wall ghost copies the
density, axial momentum, and pressure and negates the radial velocity
($v \to -v$), which makes the axis flux mass-free: by symmetry nothing
crosses the axis, and the negated ghost cancels the mass transport in
the Riemann flux. At the axis face $r_- = 0$, so the first radial row
only sees its outer face; the $r$-weighted update needs no special
treatment beyond the boundary condition.

## Validation

Mach-2 flow over a sphere: the bow shock standoff measures 0.317 against
the Ambrosio-Wortman 1962 correlation $0.143 \exp(3.24/M^2)$, which
gives 0.321 at Mach 2, a 1.2 percent miss. The post-shock density peaks
at 3.03, inside the Rankine-Hugoniot band. The planar cylinder's
standoff at the same Mach number is roughly half of this, and the
geometric source terms are what distinguish the two flows: the sphere
case is the test that they are present and correct.

The setup and the checked numbers are on
[the sphere validation page](../../verification/axisymmetric-sphere.md).

## See also

- [2D Euler solver](2d-euler.md), the planar march whose fluxes the
  annular update reuses
- [Sphere validation](../../verification/axisymmetric-sphere.md)

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
[[verification/axisymmetric-sphere|the sphere validation page]].

## The chamber lane

The axisymmetric march carries the nine-species combustion closure, so a
chamber and its throat run on the same code. Three pieces make the chamber
behave like a chamber rather than a pipe:

- The radial species transport carries the face radii, like the flow
  update does, and takes the time step on its own: the annular divergence
  already holds the cell width, so dividing by the width a second time
  moves the wrong amount of composition into the inner rows.
- The chamber's inflow is a real-gas reservoir. The ghost expands the
  reservoir state to the local face pressure through the closure, so
  the chamber holds its pressure while the throat chokes, and the
  pressure wave that the throat sends back upstream is answered
  instead of reflected. The perfect-gas reservoir ghost cannot answer
  a real-gas closure, so the combustion lane has its own. The feed
  passes through an injector, whose area the case gives, so the chamber
  takes in the mass the throat lets out.
- The chamber starts from the quasi-1D solution of its own wall
  contour, so the march begins at the state the throat will settle
  into rather than at a uniform fill it has to discover. The contour
  comes from the case's polygon, and the chain stops at the closing
  edge: the vertex that closes the polygon over the top of the
  chamber sits at the box's outer radius, and seeding the duct with
  that radius blows the exit area up by a factor of twenty.

The body fill reads the case's injection composition, not the scalar
inert reference: under the species closure the inert reference is
unburned gas at an absurd temperature, and the reconstruction stencil
would feed it to the wall-adjacent cells.

## The time step

The axisymmetric lane's radial term is stiff, and the time step has to
respect it. A Courant number of 0.05 holds a chamber together cleanly.
At 0.15 and above the radial update goes unstable: the pressure and
temperature alternate between neighbouring rows, the hot rows lose their
energy to the cold ones, and the chamber collapses to the closure's
temperature floor within a few hundred steps. The alternating pattern
looks like a geometry problem, and the face radii are fine. It is the
time step, and the cure is a smaller one.

The planar lane tolerates a larger Courant number than the annular one,
so a case that carries both kinds of flow needs the smaller of the two.

## See also

- [[numerics/2d/euler|2D Euler solver]], the planar march whose fluxes the
  annular update reuses
- [[numerics/2d/combustion|Combustion]], the nine-species closure and march
- [[verification/axisymmetric-sphere|Sphere validation]]

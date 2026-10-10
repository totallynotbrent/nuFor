---
title: Cut cells
---

The staircase mask is the cheap way to put a solid body on a Cartesian
grid: freeze every cell whose center falls inside the body, and the wall
the solver sees is the boundary of the frozen block, up to a full cell
off the true surface wherever the surface runs oblique to the grid
lines. Cut cells replace the staircase in the cut band with fractional
geometry. A cell straddling the wall carries a fluid volume fraction and
an open aperture on each of its four faces, the update sees how much of
the cell and of each face is fluid, and surface quantities sample the
true geometry instead of the mask band.

The fractions live in `crates/nufor-core/src/cutcell.rs`, the march in
`crates/nufor-core/src/cutcell_march.rs` (`advance2d_axi_cut`, built on
the axisymmetric flux machinery, so the body is a solid of revolution).
The staircase march stays untouched as the control: the cut path is a
separate step function that reads precomputed fractions, and every
measured claim below is a comparison against it.

## Fractions from the signed distance

A cut cell is described by `CutFractions`: the fluid volume fraction
`vol` and the open apertures `west`, `east`, `south`, `north`. The body
is a signed distance function, positive in the fluid, and the fractions
come from sampling it.

Two builders exist. The midpoint builder samples four points per face
and a 2x2 tensor product per cell, and it fails on thin slivers: all the
volume samples can fall inside the solid while the face samples stay
open, so the volume rounds to zero against open apertures and the update
divides by the 1e-9 clamp, a billionfold amplification on the first
step. The march uses the corner builder, `cut_fractions_corner`, which
samples the signed distance at the four cell corners, places the edge
crossings by linear interpolation, and reads the volume off the fluid
polygon (the fluid corners plus the crossings, in order, by the
shoelace area). It is exact for a planar interface and second order in
the interface position for smooth walls, and a thin sliver keeps the
positive volume its open aperture says it must have.

Two rules inside the corner builder carry the correctness:

- The edge orientation. On an edge whose endpoint samples are `a` then
  `b` of opposite sign, the crossing sits at `t = a/(a-b)` from the
  first endpoint. The fluid span is `t` when the first endpoint samples
  fluid, and `1-t` when it samples solid. Getting this backwards
  mirrors every aperture in the band.
- The zero-volume closure. A cell with no fluid volume is solid, and
  all four apertures close with it. The geometry and the update must
  agree that the cell passes nothing: a zero-volume cell holding an
  open aperture is the sliver failure again, moved from the sampling
  into the update.

Each cut cell also carries its embedded-boundary piece: the interface
length comes from the aperture imbalance, which closes the geometry the
way the discrete divergence does, and the interface normal comes from
the signed distance at the cell center, falling back to the aperture
direction when the caller provides none.

The unit contracts live in
`crates/nufor-core/tests/cutcell_geometry.rs`: the cut band is a proper
band, with no lone cut cell isolated from cut or solid neighbors; the
fractions vary smoothly along it, with no jump larger than a full cell
band; and the total solid volume matches the analytic body to better
than one cell. The builder's own tests pin the exact cases: a half-plane
cut is exact, a diagonal gives the corner triangle's exact area and
open spans, and a thin sliver keeps a positive volume with an open
aperture.

## The fractional update

One step, in divergence form:

$$
\partial_t U = -\frac{\sum_f \alpha_f A_f F_f - A_{eb} F_{eb}}{\alpha_V V} + S(r)
$$

Each open face contributes its HLLC flux scaled by its aperture, with
the annular face radii and the `p/r_c` geometric source carried over
from the axisymmetric update unchanged. The solid part of the cell
boundary contributes the embedded-boundary wall flux, and the whole
divergence sum divides by the volume fraction. That division is the
small-cell problem. The timestep is the global CFL step, the fastest
wave over the fluid cells against the full cell width, never the
fractional width. The small cells ride the stabilization instead of
shrinking the march. One step runs the shared flux sweep, forms each
fluid cell's fractional update with the wall flux joined into its
divergence sum, applies the increments, then pools the cut band
through the redistribution.

The wall flux is a reflected-ghost HLLC pair in the wall-local frame:
rotate the cell state so the axis aligns with the cell's outward wall
normal, reflect the normal velocity, solve the pair, rotate the flux
back. It enters the divergence sum with the same outward-positive sign
as the open faces, and it carries the `dt` factor: the wall impulse is
`dt` times flux times interface length, and dropping the `dt`
overscales it by `1/dt`.

The exact inviscid wall statement, zero mass, momentum `p n`, zero
energy, was tried first. It rings. A momentum-only wall update with no
dissipation oscillates, the timestep collapses toward 1e-6 while every
state stays physical, and the run dies with no NaN and no negative
density to flag it. The reflected pair's dissipation is load-bearing,
not a defect. Its cost shows up at hypersonic walls.

## The small-cell problem

The explicit update on a cell with volume fraction $\alpha_V$ is stable
only at a timestep that shrinks with $\alpha_V$. A cell holding a
thousandth of a volume at the global timestep is unstable by three
orders of magnitude, and a global march cannot afford per-cell step
sizes. Three families were built and measured on the Mach-2 sphere
(180x90, the sphere validation case):

- Improvised redistribution variants, a stashed numerator, keep-rate
  scaling, large-cell-only receivers, receiver eligibility gating: all
  failed. The excess a flux update hands off is not a gas state, it is
  mass-rich and energy-poor, it poisons the receivers' thermodynamics,
  and handoffs chained through small cells concentrate the excess
  instead of spreading it.
- Chern-Colella flux redistribution: diverged within `t = 0.7`. The
  face fluxes have already moved the gas, and the excess handoff moves
  it again, so the receivers take a double extraction.
- Berger-Giuliani state redistribution: stable through the full sphere
  run, `t = 4.0` in 17229 steps. A cut cell grows a ring of fluid
  neighbors until the pooled fluid volume reaches one full cell, each
  neighborhood takes a volume-weighted average of the updated states,
  and the final state is the average of the covering neighborhoods,
  with overlap counts so shared cells are not counted twice. The paper
  pools only cells below half a volume. In the axisymmetric setting the
  annular divergence near the axis amplifies the fractional update past
  what that half-cell rule tolerates, and the paper's threshold ran
  stable only to about `t = 1.1` on this case, so the march pools every
  cut cell.

What the stable scheme buys is the surface read. At the stagnation
point the cut march reads a density of about 4.0, the isentropic
stagnation compression, where the staircase band reads 1.9 to 3.0
because its sampled cells sit a full cell off the surface. The cost is
mixing irreversibility near the wall: the measured bow-shock standoff
shifts to 0.217 against the staircase's 0.317 on the same grid, with
the Ambrosio-Wortman reference at 0.321.

## The hypersonic limit

At Mach 22 two mechanisms defeat the scheme, both measured on the
capsule case at two resolutions:

- The wall pump. The reflected-ghost pair at a hypersonic wall is a
  mirror-state Riemann problem, not a wall reflection: its star pressure
  pushes the surface cells toward roughly twice the stagnation
  pressure. Measured on the 160-cell grid, stagnation Cp 3.91 to 3.99
  where the physical value at that resolution is about 1.7 to 1.9, and
  C_A 1.86 to 1.91 against the Newtonian 0.603. It over-reads about as
  badly as the staircase under-reads.
- The pool straddle. The redistribution neighborhoods span 9 to 15
  cells while the shock layer is 1 to 2 cells thick at practical
  resolutions, so the pooled surface reads average pre-shock and
  post-shock states.

The outcome: at 160 cells the run completes but over-reads the
surface; at 640 it destroys the field within 700 steps, with more than
ten thousand negative densities along the flank and the wake cut band.
The C_A gap does not close through this implementation, and the honest
number remains the staircase C_A of about 0.256, flat across the 320,
640, and 1280 runs, so that bias is systematic rather than a
resolution artifact. The next rung is a wall-local scheme that keeps
the shock layer out of the averaging, multi-level flux redistribution
or an h-box method, and that is research now, not a configuration
change to the current scheme.

## Conservation

On a cut mesh the conserved quantity is volume-weighted: mass is the
sum over fluid cells of `alpha * V * rho`, and energy likewise. Solid
cells hold the quiescent reference state, not gas. An unweighted cell
sum is the wrong quantity and shows a few percent of apparent drift
that is not a leak.

The contract is `crates/nufor-core/tests/cutcell_conserve.rs`: the
sphere-cone inside a closed slip-walled domain, 40 steps, with the
volume-weighted mass and energy drift asserted under 2e-2 and 3e-2. In
a closed domain the interior fluxes telescope exactly, so whatever
drift remains measures what the wall treatment injects. The band is
first-order EB wall work rather than machine precision: the
reflected-ghost wall flux does artificial work at the embedded
boundary, with a local error of order `dt / alpha_V` that the state
redistribution pools but cannot fully cancel. The real leaks met
during the build, a missing `dt` on the wall flux, a broken flux
telescope, an unbalanced redistribution, all showed up 10 to 100 times
larger.

## See also

- [[numerics/2d/axisymmetric||Axisymmetric Euler]], the annular update
  whose flux machinery the cut march reuses
- [[numerics/2d/supersonic-cylinder||Supersonic cylinder]], the
  immersed-body mask the cut band replaces
- [[formats/case-body||Immersed bodies and axisymmetric runs]], the
  `[body]` case section that puts a solid in the domain
- [[verification/axisymmetric-sphere||Sphere validation]], the Mach-2
  case the stabilization was measured on
- [[aeroshell||Aeroshell cross-check]], the Mach-22 capsule numbers

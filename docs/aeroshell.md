---
title: Aeroshell cross-check
---

The sphere-cone capsule cross-check: the user's entry-vehicle aerodynamics
(a modified-Newtonian panel method over the 10080-facet STL) against a
perfect-gas axisymmetric CFD run at the peak-q trajectory point. It is the
sprint deliverable for rungs 2 and 3 of the roadmap.

## The case

The capsule: spherical nose R_N = 1.90 m, tangent cone half-angle
9.566 deg flaring to the base radius 2.365 m at x = 4.5 m, base area
17.57 m^2. The trajectory point is peak dynamic pressure: about Mach 22
at 30.1 km, 6675 m/s. The CFD run is perfect gas (gamma 1.4), alpha 0,
on the 640x320 fine grid with the body rasterized as a solid-cell mask.
No chemistry, no heat flux: above roughly Mach 8-10 only the pressure
field is meaningful, which is exactly what this cross-check uses.

## What the CFD reproduces

- The detached bow shock ahead of the nose, with post-shock density
  peaking at 6.41 against the Rankine-Hugoniot high-Mach limit of 6.0
  (the overshoot is the MUSCL limiter at the axis).
- The stagnation Cp: 1.131 on the 320x160 grid, 1.670 on the 640x320
  grid, 2.43 on the 1280x640 grid. The 640 value converges toward the
  Rayleigh-pitot perfect-gas Cp_max of 1.838; the 1280 value
  overshoots it by a third, which pins the overshoot on the nose
  cell's placement and the limiter rather than on a physical
  mechanism (the true stagnation value cannot exceed the pitot
  value). The honest read: the surface stagnation Cp converges
  non-monotonically and grid refinement alone does not settle it.
- The user's Cp_max 1.816 (effective gamma 1.2, entry conditions)
  sits within 1.2 percent of the perfect-gas Rayleigh value: at Mach 22
  real-air effects move the stagnation pressure only slightly.

## The cross-check numbers

The modified-Newtonian axial force integrates to C_A = 0.573 with his
Cp_max over the analytic meridian, 0.580 with the perfect-gas Cp_max,
and 0.603 by his own panel integration over the STL facets (the facet
integral reproduces his table value exactly, 0.6031). The three agree
within 5 percent, which validates the panel method and the analytic
meridian against each other.

The CFD profile-integrated C_A comes out 0.257 at 640 and 0.256 at
1280 (0.247 at 320): flat across a 4x refinement, so the staircase
sampling bias is systematic, not a resolution artifact. It is
localized: 97 percent of the Newtonian axial force comes from the
outer band of the spherical cap, and the staircase mask's sampled
cells sit a full cell off the true surface exactly there, reading
the already-expanded pressure instead of the wall pressure (the sim
reads Cp near 0.19 across that band where Newtonian expects 1.1 to
1.7). The density field itself shows the shock layer hugging the
body.

Cut cells were built and tried (embedded-boundary fractions from the
signed distance field, corner-interpolation apertures, an exact
slip-wall flux, Berger-Giuliani state redistribution for the small
cells). The verdict from the mach-22 re-run: the machinery reads
true stagnation states at moderate mach (the mach-2 sphere recovers
the full isentropic stagnation compression the staircase cannot
see), but at hypersonic walls the reflected-ghost wall Riemann
problem pumps the surface cells to roughly twice the stagnation
pressure, and the redistribution neighborhoods average across a
shock layer that is one to two cells thick at practical resolutions.
At 160 cells the run completes but over-reads the surface; at 640 it
destroys the field within 700 steps. The pure pressure-force wall
flux (the exact inviscid layout) rings and collapses the timestep in
every stabilization scheme tried. So the C_A gap does not close
through this cut-cell implementation: the failure is documented with
the mechanism identified (wall pump plus pool straddle), and the
next solver rung is a wall-local redistribution that keeps the shock
layer out of the average (multi-level flux redistribution or
h-box), not a configuration change.

## Figures

- the 3D shell colored by his own panel Cp (no CFD claims on it),
- the fine-grid density field with the bow shock and the sampled
  surface,
- the Cp(x) cross-check plot with the divergence annotated.

## Honesty boundaries

Perfect gas only. No chemistry, no heat flux, no temperature claims,
alpha 0. The heating numbers for the capsule come from the user's own
trajectory tooling. Every figure states perfect gas.

## See also

- [Axisymmetric Euler](numerics/2d/axisymmetric.md)
- [Sphere validation](verification/axisymmetric-sphere.md)
- [Roadmap](roadmap.md)

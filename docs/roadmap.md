# Roadmap

nuFor's trajectory is a ladder: each rung is a validated solver
capability with its analytic gates, in an order where every milestone
makes the next one possible. The calendar is intentionally loose: the
binding constraints are the power-limited box and the session pace, not
dates. Milestones land when their gates pass.

## Where we are

- 2D Euler, viscous, and RANS (Spalart-Allmaras): validated, committed.
- Web UI: case gallery, mesh view, live runs, validation dock.
- 3D inviscid Euler: validated, committed.
- 3D viscous Navier-Stokes: validated, committed.
- Masked bodies (rung 2): done, validated on the cylinder bow shock and
  the sphere-cone trade-study geometry.
- Axisymmetric solver (rung 3): done, validated on the sphere bow shock
  against the Ambrosio-Wortman standoff. The user's capsule ran at the
  peak-q point, Mach 22 perfect gas: bow shock, post-shock states, and
  the surface Cp cross-check against his modified-Newtonian numbers.

## The ladder, in order

1. **3D viscous**: done.
2. **Masked bodies**: done. The solid-cell mask answers any shape that
   provides inside plus normal: the analytic circle, polygons, and
   analytic signed-distance bodies. Staircase walls are first-order at
   the surface but give real bow shocks, standoff, and pressure
   distributions.
3. **Axisymmetric solver**: done. The annular update carries the face
   radii and the p/r geometric source, with the slip-wall ghost as the
   axis condition.
4. **3D SA/RANS**: Spalart-Allmaras in 3D with wall distance.
   Validation: law-of-the-wall channel/duct, spanwise-symmetric plate.
   Roughly 1.5-2 weeks of sessions.
5. **Unstructured meshes**: triangles in 2D, tets in 3D. Removes the
   staircase limitation for good. The largest single lift on this list.
6. **Cut cells**: second-order accuracy at curved walls instead of the
   staircase's first order. Only worth it after unstructured work
   settles.

Open research items, worked when they block a milestone rather than on
a schedule: LTS (local time stepping) convergence, and the near-wall
u+ odd-even checkerboard from the grid-convergence study.

## Explicitly out of scope

Real-gas effects: dissociation, ionization, vibrational excitation.
nuFor is a perfect-gas solver; above roughly Mach 8-10 the pressure
field is still meaningful but temperatures and heat flux are not, and
the solver stays out of that regime. Anything hotter gets handed off to
tools built for it (LAURA-class codes), with nuFor providing the
perfect-gas cross-check. This boundary is a feature, not a gap: it is
what keeps the validation culture honest.

## Aeroshell work

The user's sphere-cone capsule work (entry trajectory, TPS trades,
structural cases) is the motivating customer for rungs 2 and 3. The
immediate sprint: masked bodies + axisymmetry + a perfect-gas
cross-check of the modified-Newtonian aerodynamics at the peak-q
trajectory point, with report renders. Real hypersonic heating numbers
for that project come from the user's own trajectory tooling, not from
nuFor.

## See also

- [Architecture](architecture.md)
- [Verification](numerics/foundations/verification.md)
- [Regression suite](verification/regression-suite.md)

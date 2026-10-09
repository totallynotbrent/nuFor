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
6. **Cut cells**: fractional-cell geometry at curved walls instead of
   the staircase, so surface pressures sample the true surface. In
   flight: the corner-interpolation fractions, the embedded-boundary
   wall flux, and Berger-Giuliani state redistribution are built and
   stable through the Mach-2 sphere run; the accuracy contract against
   Ambrosio-Wortman is the remaining gate. First order at the cut cells
   for now; the second-order flavor is a later rung.

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

## The case workflow (the product thesis)

The point of nuFor is the standard CFD loop, the way desktop tools do
it, scoped to what this solver can honestly compute:

1. **Mesh in.** Bring a geometry in (an uploaded STL watertight body,
   or an analytic body), sit it in a control volume you size in the
   case file.
2. **Set the flow.** Mark which domain faces are inflow, outflow, or
   wall; set the inflow velocity, density, pressure, and Mach regime.
   Everything is a case.toml save file: load it, edit any field, save
   it back, run it. Nothing about a case is hardcoded in the server.
3. **Run.** The march runs server-side with checkpoints; the UI shows
   live frames.
4. **See results.** Field views of pressure, speed, density,
   temperature, Mach on any plane: 2D cases whole-domain, 3D cases by
   slice (top-down, side, any xy/xz/yz offset). Surface tables: the
   pressure distribution along the body, integrated forces (the
   aeroshell's axial force coefficient comes straight off this).
   Actual flow visualization (streamlines, particles) is a later
   nicety, not a milestone.
5. **Keep it.** Save the case plus results as files, export any view
   as PNG or PDF.

The case-file system is the current milestone: the case.toml schema
gains a body section, the server gains case list/load/save/delete
routes, and the UI gains the editor plus plane-slice views. The
detailed gates live in the working goal file for the session.

## Aeroshell work, now

The user's sphere-cone capsule is the first customer of the case
workflow. The cut-cell rung (rung 6) is built and measured: at
moderate mach it reads the true stagnation states, but at the
mach-22 surface the wall Riemann pump and the redistribution pool
straddle defeat it, so the C_A gap stays open through this
implementation and the honest number remains the staircase 0.256.
The next solver rung toward the surface is a wall-local
redistribution (multi-level flux redistribution or h-box) that
keeps the shock layer out of the averaging. Real hypersonic heating
numbers for that project come from the user's own trajectory
tooling, not from nuFor.

## Rocket engine CFD, later

The eventual customer the case workflow is being built for. When the
aeroshell scope closes: nozzle and internal flows over uploaded
geometry, still perfect-gas. Combustion, chemistry, and heat transfer
stay beyond the solver's honest envelope; the engine work is the
cold-flow aerodynamics and the perfect-gas cross-check, the same role
nuFor plays for the aeroshell.

## See also

- [Architecture](architecture.md)
- [Verification](numerics/foundations/verification.md)
- [Regression suite](verification/regression-suite.md)

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
  against the Ambrosio-Wortman standoff. The capsule ran at the
  peak-q point, Mach 22 perfect gas: bow shock, post-shock states, and
  the surface Cp cross-check against the modified-Newtonian numbers.
- Equilibrium air (CEA closure): done. The planar and axisymmetric 2D
  paths run the CEA table under `eos = "eqair"`; the legacy TGAS1 fits
  remain as the face-state extrapolation the wake needs. A four-station
  trajectory sweep (30/34.6/40/42.5 km, alpha 22.5) is validated at
  320x176 and anchored at 1280x640.
- Performance: the static body classification is cached and the closure
  passes are threaded; the hi-res anchor runs 8x faster than the serial
  march (3.96 h for 11,549 steps at 1280x640 on the 35 W box).
- Fay-Riddell stagnation heating gate (rung 6): done. The
  wall-resolved cylinder gate closes at conductor 10.7 W/cm2 vs the
  planar closed form 8.5 on the march's own edge state (+27%, the
  conductor still converging downward); wall heat flux is citable
  output now.
- Unstructured 2D: a gmsh-driven finite-volume path exists
  (`ugrid.rs`, advance_ugrid) with the Sod gate; a prototype, not yet a
  validated rung.

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
4. **Equilibrium air**: done (moved up from "out of scope": the
   aeroshell project needed it, and the CEA table plus the TGAS1 face
   fits are the validated form).
5. **Live two-temperature nonequilibrium (Park 2T)**: done. `e_v` is a
   live field (unpinned t_v) with Millikan-White + Park-cutoff
   relaxation and the t_rxn = sqrt(T Tv) mixing rule; the 5-species
   rates (Park 1990) sit in park.rs waiting for rung 8. Gates: the
   analytic Landau-Teller decay box, energy conservation to machine
   precision, then the Mach 20 nitrogen cylinder (Casseau et al. 2016)
   as the external anchor. License note: hyStrath is GPL-3.0, same as
   nuFor; the borrow is the validation target and the physics
   constants, not code.
6. **Fay-Riddell stagnation heating gate**: done. A
   wall-resolved viscous cylinder at entry conditions, validated
   against the planar Fay-Riddell closed form (0.57 Pr^-0.6 vs the
   sphere's 0.763). The gate case runs a nose-only domain at global
   dt: local time stepping detonates on this class at long times
   (documented on the viscous page), and the nose-only cut holds the
   march cost to hours on this box. The wall closures it exercised:
   the reflected-riemann slip wall in the advective stage, the
   one-sided conductor in the viscous stage, and closure-lane
   conduction with real units. The planar closed form is in the tree
   (`fay_riddell_q_2d`), pinned to the Beckwith-Cohen similarity
   constant and the Rankine-Hugoniot tables. This is the gate that
   makes any wall heat flux the solver prints citable.
7. **Continuum-validity probe per station**: Knudsen number and Bird's
   breakdown parameter from the freestream, one line in the run log.
   Cheap, and it puts a number behind "continuum CFD applies here"
   instead of an assertion.
8. **Species continuity (5-species Park nonequilibrium)**: the big
   real-gas rung. Species mass fractions as transported fields with the
   Park 1990 dissociation/exchange set already tabulated in park.rs,
   D0 energies routed through the closure, positivity handling, and a
   nonreacting path that pays no reacting overhead. Gates: mass
   conservation to roundoff, equilibrium composition against published
   air tables at fixed (T, p), then the Mach 20 nitrogen cylinder with
   dissociation. The reacting-flow ordering in the main README holds:
   this rung starts only after thermodynamics and transport are
   separately trustworthy.
9. **Catalytic wall brackets**: non-catalytic and fully-catalytic wall
   options for the diffusive species flux at the surface. With rung 8
   this brackets the heating answer for TPS (fully-catalytic is the
   conservative bound a report wants). The Maxwell velocity-slip and
   Smoluchowski temperature-jump wall conditions from the same body of
   work are optional extras: below 43 km the trajectory is continuum
   (Kn << 1), so slip walls are a research branch, not a rung.
10. **3D SA/RANS**: Spalart-Allmaras in 3D with wall distance.
    Validation: law-of-the-wall channel/duct, spanwise-symmetric plate.
    Stays low in the list: the entry trajectory is laminar-dominant,
    and hyStrath's own turbulence story is unchanged OpenFOAM defaults,
    so there is nothing to borrow.
11. **Adaptive mesh refinement**: block-structured AMR on a density or
    Mach-gradient indicator, in the hyStrath DyM spirit. The hi-res
    anchor spends four uniform-grid hours on a shock layer occupying a
    few percent of the cells; two AMR levels put that under an hour on
    this box. Placed after 2T/species because the indicator should
    follow the nonequilibrium structure, not just the shock.
12. **Unstructured meshes**: triangles in 2D, tets in 3D. The gmsh
    prototype exists and the research note is written; promoting it to
    a validated rung is the largest single lift. Arbitrary geometries
    without the staircase are the payoff: the structural advantage
    OpenFOAM-based solvers inherit for free.
13. **Cut cells**: fractional-cell geometry at curved walls instead of
    the staircase, so surface pressures sample the true surface. Built
    and measured: true stagnation states at moderate Mach, but at
    hypersonic walls the wall Riemann pump and the redistribution pool
    straddle defeat it. The next solver step here is a wall-local
    redistribution (multi-level flux redistribution or h-box) that
    keeps the shock layer out of the averaging.
14. **MPI-parallel Fortran core**: the parallelism rung. Fortran owns
    the parallel solver core (domain decomposition, halo
    exchange, the march) behind one coarse FFI; Rust stays the front
    door (case files, CLI, UI, output) and never calls MPI. Design,
    gates, and measured scaling: [the MPI page](numerics/parallel/mpi.md).
    Phase F1 parallelizes the euler_2d march (the capsule path) with the
    np=1 bit-identity anchor; phase F2 adds viscous+SA; phase F3 lands
    all new physics (live 2T, species) MPI-native from day one, which
    reorders rungs 5 and 8 to sit inside this one.

## What we deliberately do not borrow

- MHD (Lorentz/Joule source terms, conductivity models, Hall effect):
  the aeroshell trajectory has no magnetic interaction. This is
  hyStrath's differentiator and dead weight here.
- DSMC and hybrid CFD-DSMC: below 43 km the flow is continuum by orders
  of Knudsen. The breakdown probe (rung 7) is the honest boundary.
- 11-species ionized air and electron-energy modes: ionization is
  negligible below ~9.5 km/s at these altitudes; five species covers
  air chemistry to about 12,000 K.
- Two-equation turbulence models (k-epsilon, SST): SA exists; the
  laminar entry regime does not need a model zoo.
- The OpenFOAM framework itself: GPL-compatible but a build and runtime
  dependency from another era (OF v1706); nuFor's value is being
  self-contained, inspectable, and honest about its envelope.

Open research items, worked when they block a milestone rather than on
a schedule: LTS (local time stepping) convergence, and the near-wall
u+ odd-even checkerboard from the grid-convergence study.

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

The four-station trajectory sweep at alpha 22.5 is validated at
320x176. The honest reads: L/D 0.55-0.72 rising with altitude across
the sweep; peak shock-layer temperatures 5.1-14.0 kK confirm the
thermal environment scale; peak pressures to 1.09 MPa now read off the
rebuilt CEA grid (101 pressures to 10.1 MPa, 30 temperatures to
30000 K, sensible internal energy) instead of clamping at the old
table's 1 MPa cap. The 1280x640 anchor re-measured on the new grid
reads C_D 0.406, L/D 0.371, peak pressure 1.17 MPa at 8217 K, and
lands inside the Newtonian trim band; the previous anchor numbers
(3.96 h for 11549 steps) predate the rebuild. Wall heat flux is
not citable until the Fay-Riddell gate (rung 6) passes. The defensible
report posture: Newtonian numbers stay the design basis, the CFD sweep
is the confirmatory run with a limitations paragraph, and rungs 5-6
are the upgrade path.

## Rocket engine CFD, later

The eventual customer the case workflow is being built for. When the
aeroshell scope closes: nozzle and internal flows over uploaded
geometry. Combustion chemistry arrives only after rung 8's species
continuity is validated; until then the engine work is the cold-flow
aerodynamics and the perfect-gas cross-check, the same role nuFor
plays for the aeroshell.

## See also

- [Architecture](architecture.md)
- [Verification](verification/index.md)
- [Regression suite](verification/regression-suite.md)
- [Equilibrium air](numerics/foundations/equilibrium-air.md)

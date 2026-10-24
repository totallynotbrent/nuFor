---
title: 2D boundary conditions
---

The 2D solver exposes a boundary condition on each of the four sides. A ghost
cell is filled per side before the MUSCL reconstruction, so the flux step is
unchanged.

## Conditions

- Transmissive: the ghost copies the interior cell, and waves leave freely.
  This is the default and the only behaviour the solver had before this step.
- SlipWall: the ghost mirrors the velocity component normal to the wall and
  copies density, pressure, and the tangential velocity. A closed slip box
  conserves mass exactly and a uniform stream passes over a slip wall
  undisturbed (both asserted in the tests).
- NoSlipWall: both velocity components mirror, so the wall-average velocity
  vanishes. The viscous pass adds the one-sided wall shear, so a boundary
  layer develops against the wall.
- SupersonicInflow: the ghost is held at a fixed freestream primitive state,
  which a supersonic face cannot send a signal back into.
- SupersonicOutflow: the ghost copies the interior, the standard exit for a
  supersonic face.
- ProfileInflow: the ghost carries a boundary-layer profile instead of one
  uniform state. The profile prescribes u(y) and v(y) along the inflow plane
  (evaluated at each row's true height, so it is correct on clustered meshes),
  while density and pressure stay at the freestream values. Two sources: the
  built-in Blasius laminar layer anchored at a virtual leading edge, or a
  tabulated y/u/v file. The viscous pass uses the same profile for its ghost
  gradient, so the incoming shear the layer feels matches what it carries.
- ReservoirInflow: a subsonic stagnation-state inflow. The ghost expands
  (pt, Tt) isentropically to the local interior face pressure, so the inflow
  is not a fixed state: a compression wave traveling upstream meets a ghost
  that answers with a denser, slower state, the acoustic feedback a
  supersonic-inflow ghost cannot give. The tangential velocity is zero.
  This is the chamber-inflow condition the nozzle case runs on.
- ReservoirInflowRealGas: the same idea for a real-gas lane, where the
  stagnation state is a composition rather than a temperature. The ghost
  takes the interior cell's own pressure and temperature, the state the
  face sees, and the velocity that follows from the reservoir's expansion
  to that pressure. The mixture's own gamma comes from the closure at the
  reservoir temperature, not from the shaping constant. A combustion case
  selects it by asking for a `reservoir_inflow` on the left side; the
  composition and the reservoir state come from the case's injection block.

  The ghost's temperature is the chamber's, not the reservoir's. The gas
  leaving a reservoir has expanded and mixed into the chamber by the time
  it reaches the face, and the reservoir's own expansion runs colder than
  the chamber's gas. That difference is not a rounding error: the
  combustion closure's internal energy carries the formation offset, so a
  temperature difference between the ghost and the cell is a megajoule
  scale enthalpy difference. The face's energy flux would carry that
  difference into the first column, the column would cool, its pressure
  would fall, the ghost would thin, and the feed would run away.
- An injector area: a chamber is fed through an injector at its head, not
  across its whole face, so the reservoir inflow takes an optional
  `injector_area` in square metres. Set it in the case:

        [boundaries]
        left = "reservoir_inflow"
        injector_area = 7.2e-4

  An injector is a converging nozzle, and that shapes the whole feed. The
  mass flux through it follows the nozzle relation between the tank state
  and the chamber's pressure. Once the chamber sits below the tank by more
  than the critical pressure ratio, the flow reaches a sonic exit and
  chokes: the mass flow then follows the tank state and the injector's
  area alone, and the chamber's pressure can no longer slow it down. Below
  the critical ratio the injector runs unchoked and its flow follows the
  chamber's pressure.

  The feed's mass flow is the injector's area times the mass flux, spread
  over the face the injector enters, so the chamber balances the throat
  instead of either stalling or flooding. Nothing about the geometry is
  assumed: the critical ratio and the choked flux come from the mixture's
  own gamma and the tank's own state.

  The injector sets the mass flow and leaves the ghost's state alone. A
  ghost carrying the injector's own exit state would sit at the sound
  speed against a slow chamber, and the cells next to the face would pile
  up mass and drain their energy trying to match it. The injector's flow
  reaches the chamber through the face, and the ghost's velocity is what
  carries that flow.

  The face's flow area counts the fluid rows only. A chamber's wall
  covers the rows above it, and those rows carry no flow, so counting
  them would spread the injector's flow over a face the gas never sees
  and quietly starve the chamber.

  The tank's density has to match its pressure and temperature. A tank at
  120 bar and 3677 kelvin with a mixture constant of 527 j per kilogram
  kelvin holds about 6.19 kilograms per cubic metre, and a tank density
  that belongs to a different pressure quietly starves the chamber.

  A chamber that runs on this holds its pressure. On a 240 by 60 grid
  with a chamber, a throat and a bell, marched 3000 steps at a Courant
  number of 0.05, the chamber sits at 66.1 bar against the 70 bar its
  seed carried, the throat sits at 40.5 bar against the 40.54 bar an
  equilibrium calculation gives for the same mixture, and no cell in
  the domain falls below 2000 kelvin. That run's injector was 7.2
  square millimetres, which passes the 4.01 kilograms per second a
  choked throat of that mixture swallows.
- SubsonicOutflow: the ghost copies the interior but pins the ghost
  pressure to a fixed back value, the standard subsonic exit for the
  pressure-ratio nozzle studies. The config surface is
  `[boundaries] right = "subsonic_outflow"` with the back pressure read
  from `[boundaries.outflow_state]`.

The two engine conditions live on the rust lanes (the axisymmetric and
planar marches); the MPI Fortran core carries them as transmissive until
its boundary vocabulary grows.

## Use

    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow { rho: 1.0, u: 1.5, v: 0.0, p: 1.0 },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };
    advance2d_rk2(&mut state, &grid, gamma, cfl, muscl, &bc)?;

A profile inflow is built once and shared by reference, since the condition
carries the profile data:

    let profile = BlasiusProfile::new(u_inf, nu, x_le)?.anchored_at(x_in);
    let bc = Boundaries2d {
        west: Bc2d::ProfileInflow {
            profile: Arc::new(InflowProfile::Blasius(profile)),
            rho: 1.0,
            p: 1.0,
        },
        east: Bc2d::Transmissive,
        south: Bc2d::NoSlipWall,
        north: Bc2d::Transmissive,
    };

## What the profile condition is not

ProfileInflow is an inflow-plane prescription, not a far-field Riemann
condition. It assumes the flow entering the domain is known, which is true
for the laminar flat plate and for profiles measured or computed upstream.
It is not a place to enforce a desired solution the interior disagrees with;
the corner where the profile meets a wall or an outflow side always carries a
small compatibility transient, and validation assertions should treat those
stations accordingly.

## Verification

- A closed slip box keeps every cell exactly at its initial state and leaks
  no mass.
- A supersonic stream through a slip-wall channel with inflow/outflow stays
  uniform to machine precision, so the walls introduce no spurious waves.
- The Blasius hold: the analytic layer laid into a clustered grid and fed by
  its own profile at the inflow is held to about two percent, and the
  discrete skin friction matches the 0.664/sqrt(Re_x) correlation within a
  few percent over the interior stations. See
  [[numerics/2d/blasius|the Blasius validation page]].

## See also

- [[numerics/2d/euler|2D Euler solver]]
- [[numerics/2d/blasius|Blasius flat plate]]
- [[numerics/2d/viscous|Viscous terms]]
- [[numerics/2d/oblique-shock|Oblique shock]]
- [[numerics/2d/shock-reflection|Shock reflection]]
- [[numerics/2d/supersonic-cylinder|Supersonic cylinder]]

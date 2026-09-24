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
  [[numerics/2d/blasius||the Blasius validation page]].

## See also

- [[numerics/2d/2d-euler||2D Euler solver]]
- [[numerics/2d/blasius||Blasius flat plate]]
- [[numerics/2d/viscous||Viscous terms]]
- [[numerics/2d/oblique-shock||Oblique shock]]
- [[numerics/2d/shock-reflection||Shock reflection]]
- [[numerics/2d/supersonic-cylinder||Supersonic cylinder]]

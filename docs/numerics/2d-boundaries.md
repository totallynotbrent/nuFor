# 2D boundary conditions

The 2D solver exposes a boundary condition on each of the four sides. A ghost
cell is filled per side before the MUSCL reconstruction, so the flux step is
unchanged.

## Conditions

- **Transmissive** — the ghost copies the interior cell; waves leave freely.
  This is the default and the only behaviour the solver had before this step.
- **SlipWall** — the ghost mirrors the velocity component normal to the wall
  and copies density, pressure, and the tangential velocity. A closed slip box
  conserves mass exactly and a uniform stream passes over a slip wall
  undisturbed (both asserted in the tests).
- **SupersonicInflow** — the ghost is held at a fixed freestream primitive
  state, which a supersonic face cannot send a signal back into.
- **SupersonicOutflow** — the ghost copies the interior, the standard exit for
  a supersonic face.

## Use

    let bc = Boundaries2d {
        west: Bc2d::SupersonicInflow { rho: 1.0, u: 1.5, v: 0.0, p: 1.0 },
        east: Bc2d::SupersonicOutflow,
        south: Bc2d::SlipWall,
        north: Bc2d::SlipWall,
    };
    advance2d_rk2(&mut state, &grid, gamma, cfl, muscl, &bc)?;

## Verification

- A closed slip box keeps every cell exactly at its initial state and leaks no
  mass.
- A supersonic stream through a slip-wall channel with inflow/outflow stays
  uniform to machine precision, so the walls introduce no spurious waves.

## See also

- [[numerics/2d-euler|2D Euler solver]]
- [[numerics/oblique-shock|Oblique shock]]
- [[numerics/shock-reflection|Shock reflection]]

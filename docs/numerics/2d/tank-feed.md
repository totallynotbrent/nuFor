---
title: Propellant tank blowdown
---

An engine does not start against an infinite reservoir. The tank is a
finite vessel of pressurised vapour, it empties as the injector drinks
from it, and the pressure the injector sees falls with it. The chamber's
answer changes as that happens, because the driving pressure for the
whole engine cycle disappears while the engine is still running.

This page records the blowdown coupling: what the tank law is, how it is
wired to the march, and how it is checked.

## The law

The tank is a rigid vessel of volume `V` holding vapour of density `rho`
at temperature `T`, discharging through an injector of area `A` against
a downstream pressure `p_down`. Its state advances in two equations.

The mass flux depends on which regime the injector is in. While the
chamber sits below the tank's pressure the injector is choked, and the
flux follows the tank state alone:

```
flux = rho * sqrt(gamma * R * T) * (2 / (gamma + 1))^((gamma + 1) / (2 (gamma - 1)))
```

Once the chamber catches up with the tank the expansion ratio falls below
critical and the injector becomes unchoked, with the flux following the
standard orifice relation at the exit plane. Between the two, the
transition is smooth.

The tank's own state follows the blowdown:

```
dm    = flux * A * dt
rho'  = rho - dm / V
T'    = T * (rho'/rho)^(n - 1)
```

where `n` is the polytropic index: 1 is isothermal, `gamma` is
adiabatic. The temperature relation makes `p / rho^n` a constant of the
motion, and the gate checks that invariant on every recorded step.

Mass drawn from the tank reduces the tank's level by exactly the same
amount, so the drawn mass and the level drop agree by construction unless
something at the seam manufactures mass.

## Wiring

The `physics.tank_feed` case stanza switches the coupling on:

```toml
[physics.tank_feed]
volume = 0.00025
rho0 = 120.0
t0 = 300.0
polytrope = 1.15
r_gas = 518.3
```

The injector's area comes from the case's `boundaries.injector_area`,
shared with the chamber's own mass-flow treatment so the two seams
compute the same number.

On every march step the solver averages the chamber-side pressure over
the inflow column, steps the tank against it for the elapsed time, and
mirrors the resulting stagnation state back into the reservoir ghost on
the west face. The chamber therefore sees a decaying reservoir rather
than a constant one. A case that omits the stanza keeps the historic
constant reservoir.

The marched case is `cases/tank-feed`: the chamber-throat geometry and
seed with the tank attached, marched at 40x16 to 1530 steps.

## What the marched case shows

Over the run the reservoir fell from 1.86588e7 Pa to 1.32552e7 Pa, a
29.0 percent drop, and delivered 7.716e-3 kg of the 3.0e-2 kg fill,
25.7 percent of the tank, in 322 seconds of wall time.

The injector ran choked on the first 428 steps and unchoked for the
remaining 1102. The chamber's pressure rose over the march, and once it
approached the falling tank pressure the injector stopped choking. The
law handles both regimes, and the gate adjudicates each step against the
regime the recorded pressures imply rather than assuming one.

| quantity | value |
|---|---|
| marched steps | 1530 |
| initial tank pressure | 1.86588e7 Pa |
| final tank pressure | 1.32552e7 Pa |
| initial tank density | 120.0 kg/m3 |
| final tank density | 89.1358 kg/m3 |
| mass drawn | 7.716061e-3 kg |
| reservoir fill | 3.0e-2 kg |
| fraction drained | 25.7 % |
| choked steps | 428 of 1530 |
| unchoked steps | 1102 of 1530 |
| achieved flux against the replay | 1.8e-13 relative, worst step |
| polytropic invariant error | 2.3e-15 relative |

## The gate

`cases/tank-feed/gates.py` replays the law independently, outside the
solver, driven by the recorded chamber pressure and timestep sequence.
It compares the replay against the solver's own record row by row:

- the flux the solver achieved on each step matches the flux the
  governing regime prescribes, to roundoff;
- the recorded tank pressure and density match the replay on every row;
- the drawn mass equals the tank's own level drop;
- the polytropic invariant holds along the record;
- the reservoir actually drained.

The gate reads the committed record in `cases/tank-feed/record/`, so it
runs in continuous integration without the five-minute march. Marching
the case locally refreshes the record and the gate prefers the fresh one.

A tank step that fails ends the run. The margin for an error there was
previously swallowed, which left the reservoir pumping mass no reservoir
could supply. The failure surfaces as a nonzero exit status and the
printed line naming the time the step broke.

## See also

- [[numerics/2d/chamber-throat|Chamber to throat]]
- [[numerics/2d/combustion|H2/O2 combustion]]
- [[numerics/2d/film-injection|Blown-in film injection]]
- [[verification/index|Verification index]]

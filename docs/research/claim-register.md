---
title: Claim register
---

Every quantitative claim the published docs make, and where it comes from.
Each row names the artifact that produces the number, so a reader can tell
whether a value was measured, gated, or taken from a reference. Values marked
provisional are not publishable as results.

## Aeroshell station sweep

Published in `research/results.md` and the equilibrium-air page.

| claim | value | source | state |
|---|---|---|---|
| station count | four, 30.0 to 42.5 km | `cases/aeroshell-h*-a22.5` | measured |
| angle of attack | 22.5 degrees | case files | measured |
| grid | 320x176 | case files | measured |
| march end time | 15 ms | run logs | measured |
| wall temperature | 300 K | case files | measured |
| grids tested against the wall-treatment defect | 120x66 and the native 320x176 | the station re-march | measured |
| density ratio, peak p, max T per station | 11.26 / 3.58e5 / 5079 and peers | `cases/aeroshell-station/gates.py` | gated |
| lift-to-drag per station | 0.372 / 0.391 / 0.381 / 0.386 | run logs | gated |
| peak wall flux per station | 0.0150 / 0.0260 / 0.0441 / 0.0766 MW/m^2 | run logs, shared estimator | gated |
| where the applied peak sits | on the hottest wall-adjacent cell of the field, 5014 K against a field maximum of 5079 K | the marched field | measured |
| peak gradient against the reported flux | consistent with the air conductivity over one cell | the marched field | measured |
| station steps and wall time | 1398 to 2647 steps, 113 to 271 s | run logs | measured |

The gate re-measures the density ratio, peak pressure, maximum temperature and
wall flux off the marched field and holds them against `reference.json` within
12 %. The reference is regenerated from a fresh sweep, never hand-edited.

## Trim anchor

Published in `research/results.md`.

| claim | value | source | state |
|---|---|---|---|
| grid and wall | 384x192, 1600 K radiation-equilibrium wall | `cases/aeroshell-capsule-trim` | measured |
| march length | 1804 steps to t = 2 ms | run log | measured |
| peak surface pressure | 0.959 MPa at (1.549, 1.281) | run log | measured |
| maximum field temperature | 25670 K | run log | measured |
| peak wall heat flux | 0.575 MW/m2 | run log, shared estimator | measured |
| C_D and L/D | 0.501 and 0.335 | run log | measured |

## Catalytic wall accommodation

Published in `research/results.md` and the catalytic-walls page.

| claim | value | source | state |
|---|---|---|---|
| non-catalytic peak wall flux | withheld | `cases/cylinder-m20-park5-visc/case-none.toml` | withheld |
| fully-catalytic peak wall flux | withheld | `cases/cylinder-m20-park5-visc/case-full.toml` | withheld |
| catalytic amplification | withheld | the two runs above | withheld |
| where the applied peak sits | the face at (0.963, 0.315) | the run log's peak line | measured |
| that cell's temperature | 18342 K, next to the field maximum of 18791 K | the marched field | measured |
| first-row temperature at the nose | 6678 K | the marched field | measured |
| first-row temperature at the shoulder | 8586 K | the marched field | measured |
| first-row temperature at the lee | 11233 K | the marched field | measured |
| temperatures of the cells beside the peak | 16137 K above, 13180 K beside, 18791 K below | the marched field | measured |
| the two implementations on the same face | 0.240 MW/m2 from the outside reporter and 3.75 MW/m2 from the march, at (0.963, 0.315) | two runs of the same case on different far fields | not comparable |
| the march's peak on the two wall conditions | 3.75 non-catalytic and 4.89 fully catalytic MW/m2 | the two runs' run logs | withheld |
| the accommodation series | 3.746, 3.820, 3.895, 4.230 and 4.893 MW/m2 at accommodation 0, 0.25, 0.5, 0.75 and 1 | the five sweep run logs | model study |
| the amplification across that ladder | 1.306 at full accommodation | the same five run logs | model study |
| the endpoints against the shipped pair | reproduced exactly | the gate in the case directory | gated |
| the monotonicity across that ladder | the peak grows with the accommodation | the gate in the case directory | gated |
| the conductance the march uses | the Wilke mixture of the RP-1232 species conductivities, varying with temperature and composition | `crates/nufor-core/src/species_props.rs` | measured |
| the earlier constant-conductance claim | withdrawn, the species lane has no constant conductance | this review | withdrawn |
| continuum validity near the wall | gradient length Knudsen number up to 1089, about 10 percent of cells above the threshold | the continuum probe | measured |
| the earlier boundary-cell attribution | withdrawn, the hot region is smooth over several cells | this review | withdrawn |
| cause of the withheld peak | the flow at the wall is out of the continuum range, so the wall flux is not a continuum quantity there | the continuum probe | diagnosed |

The earlier 0.240 and 1.216 MW/m^2 pair came from a reporter that recomputed
the physics outside the march, with a closure temperature the march never
used, on a stream the case never declared. The 0.240 figure is not
comparable with the march's own peak, and the claim that tied the two to a
constant 4.93 W/m.K conductance is withdrawn: the species lane conducts
with the temperature dependent mixture conductivity. The peak stays
withheld until the case clears the continuum criterion, and the
accommodation coefficient is now the input a measured surface
efficiency would supply.

## Stagnation heating gate

Published in `research/results.md` and the Fay-Riddell page.

| claim | value | source | state |
|---|---|---|---|
| edge plateau rho, T, p | 2.3e-4 kg/m3, 5040 K, 505 Pa | `cases/cylinder-stagnation-fay` | gated |
| converged conductor band | 10.7 to 11.2 W/cm2 | the refined 132x98 mesh | gated |
| spread of the settled band | 0.5 W/cm2 | the same mesh | gated |
| corrected closed form | 10.69 W/cm2 | CEA mixture transport at the edge state | computed |
| agreement | within 5%, and the closed form sits at the lower edge of the measured band | band against closed form | computed |

`cases/cylinder-stagnation-fay/gates.py` recomputes the edge plateau, the
settled band of the 660k-step chain, the narrower band of the refined
mesh and the ratio of the settled mean to the closed form, and
cross-checks the committed series against the run's own chunk table.

## Dissociating nitrogen cylinder

Published in `research/results.md` and the five-species page.

| claim | value | source | state |
|---|---|---|---|
| shock standoff | 0.30 m | `cases/cylinder-m20-park5` | gated |
| stagnation-line peak T | 15.3 kK | same | gated |
| stagnation-line N mole fraction | 0.34 | same | gated |
| peak N mass fraction | 0.50 | same | gated |
| drag from the pressure integral | not citable | same | provisional |
| far field the run flies | 220.0 K, 0.890 Pa, Mach 20.00 | the gate in the case directory | gated |
| far field the case declares | 220 K, 0.89 Pa, Mach 20 | `cases/cylinder-m20-park5/case.toml` | declared |

Surface-integrated drag picks up post-shock pressures where the wrapped shock
crosses the surface cells. The doc states the number is measurement-limited
rather than citable, and that stands.

The family flies the state the case declares, and the gate says so: the
far field must sit within 2% of the declared temperature and the
shock-layer scalars within a few percent of the last run. Getting there
took two fixes. The uniform seed carried the air mixture's energy into a
pure-nitrogen march, and the supersonic inflow handed the boundary pressure
to the lane's internal-energy slot, which pinned the far field to the 50 K
floor. The whole family ran on a stream nobody asked for, and the enthalpy
survived to within 0.6%, so the shock-layer scalars sat near their
published bands anyway. With the far field repaired, the standoff and the
peak nitrogen mass fraction hardly moved, the peak temperature rose 0.8
kK, and the wall-side atomic nitrogen eased to 0.34.

The gate reads `samples.csv`, the slices `samples.py` distils out of the
marched field, so the claim can be checked without marching the case
again.

## Nozzle cold flow

Published in `research/results.md` and the nozzle-cold-flow page.

| claim | value | source | state |
|---|---|---|---|
| pe/pt = 0.89 exit p and Mach | 0.874 vs 0.890, 0.372 vs 0.411 | `cases/nozzle-cdv` | gated |
| pe/pt = 0.16 exit Mach and mass flux spread | 1.742 mass-average vs 1.854, 0.07% | same | gated |
| pe/pt = 0.75 shock stand-off | shock at x ~ 6 against 1D 7.56 | same | gated |

## Scale and momentum

Not published, tracked here so the distinction is explicit.

| quantity | note |
|---|---|
| axial skin friction, wave drag | not modeled; lift-to-drag is the axisymmetric result only |
| free-electron precursor | not modeled |
| compressible base pressure | the cavity plume fills the wake, so base pressure is the freestream value |
| centrifugal correction | none; axisymmetric without swirl |
| induced drag | not computed |

## Equilibrium-air closures

Cells read the CEA table, the legacy face fits read the older surface fit.
Both are live in the same marches, so they are held to a shared tolerance.

| claim | value | source | state |
|---|---|---|---|
| agreement over the campaign envelope | median about 5 %, ninetieth percentile about 10 % | `crates/nufor-core/src/closure_tests.rs` | gated |
| envelope | 1.4e5 to 3e7 J/kg of specific internal energy, 1e-4 to 1 kg/m3 | the same gate | measured |
| hot-corner divergence above 3e7 J/kg | median about 24 % | the same gate | gated |
| cold states below 1.4e5 J/kg | outside the shared domain: the CEA table floors at 200 K while the legacy fit extrapolates | the table's own limits | measured |
| pressure to energy inversion through the legacy closure | misses by up to 100 % at the ninetieth percentile of the envelope | measurement against the campaign lattice | measured |

The gate keeps the two closures tied to each other over the region the
campaigns occupy and bounds the hot-corner departure, so a change in either
closure turns the gate red instead of quietly moving a published number. The
last row is why the legacy inversion is confined to case initialization: the
march reads the closure, it does not invert it.

## How a claim earns its state

A value is measured when it is read off a marched field or a run log, and
gated when a shipped script re-measures it and holds it against a stored
reference within a documented tolerance. Provisional values are drafts held
out of publication until a gate backs them. Numbers copied from an external
paper are labelled published reference and are not claimed as results.

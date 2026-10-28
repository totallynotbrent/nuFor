---
title: Research results
---

The numbers behind the solver features, measured on the cases I ship. Every
figure here regenerates from the case run logs in the repo, so the numbers age
with the code instead of drifting from it.

## Aeroshell station sweep

Four trajectory stations at 22.5 degrees angle of attack, equilibrium air,
viscous with a 300 K wall, 320x176 grid, marched to t = 15 ms. This is the
sweep that anchors the hypersonic 2D path, from the 30 km high-q regime
up to the 42.5 km edge of the corridor. March steps and wall time per
station are
plotted on [[numerics/foundations/equilibrium-air|the equilibrium-air page]].

| station | steps | wall (s) | max density ratio | peak p (Pa) | max T (K) | peak wall flux (MW/m^2) |
|---------|-------|----------|-------------------|-------------|-----------|---------------------------|
| 30.0 km | 1398  | 113.2    | 11.26             | 3.58e5      | 5079      | 0.0150                  |
| 34.6 km | 1999  | 155.2    | 15.79             | 4.81e5      | 7191      | 0.0260                  |
| 40.0 km | 2488  | 223.6    | 20.62             | 4.52e5      | 11611     | 0.0441                  |
| 42.5 km | 2647  | 270.8    | 18.94             | 3.68e5      | 14710     | 0.0766                  |

The march length now tracks the shock layer rather than a table ceiling.
The rebuilt equilibrium-air grid (101 pressures to 10.1 MPa, 30 temperatures
to 30000 K, sensible internal energy) carries every station without clamping,
and the shock-layer temperatures rise with altitude as the thinner air slows
the deceleration. The 42.5 km station pays the most: its thin, hot shock
layer needs the smallest time steps, and its 14700 K peak temperature is the
sternest test of the march's positivity.

Lift-to-drag across the sweep runs 0.372, 0.391, 0.381, 0.386 from the
30 km station to the 42.5 km edge: the shell holds its attitude ratio as the
freestream thins, and the wall heat flux climbs from 15 to 77 kW/m2 as the
boundary layer heats up.

The wall heat flux is the march's own number. The viscous pass evaluates the
isothermal wall conduction on the same one-sided quadratic it applies to the
field, and the run log and the progress heartbeat both read it from that one
implementation, so the two cannot disagree and no post-hoc recomputation
invents a second value. The quantities above come from the shipped gate in
`cases/aeroshell-station`, which re-measures the density ratio, peak
pressure, maximum temperature and wall flux off the marched field and holds
them against `reference.json` within 12 %.

## Stagnation heating gate

The Fay-Riddell case (`cases/cylinder-stagnation-fay`): the Casseau
Mach-20 nitrogen freestream on a 1 m cylinder nose, 1000 K isothermal
wall, nose-only domain at global dt.

| quantity | value |
|---|---|
| edge plateau $\rho_e$, $T_e$, $p_e$ | 2.3e-4 kg/m3, 5040 K, 505 Pa |
| wall conductor, converged (refined 132x98 mesh) | 10.7 to 11.2 W/cm2 |
| spread of the settled band | 0.5 W/cm2 across the band |
| closed form, pure-N2 Sutherland transport | 8.48 W/cm2 |
| closed form, CEA mixture transport at the edge state | 10.69 W/cm2 |
| agreement, converged march vs corrected closed form | within 5%, and the closed form falls inside the measured band |

The 660k-step chain on the original 84x384 grid never settled: the
conductor held a 7.4-14.0 W/cm2 band indefinitely. The tangential
direction was the problem: the Hiemenz layer is ~11 mm and the
uniform 31 mm spacing aliased it to a third of a cell. The refined
132x98 grid (100 um wall cells, 2 mm at the stagnation line) kills
the breathing entirely: the conductor settles into a 10.7-11.2 W/cm2 band
half a watt per square centimetre wide and stays there across the settled
window.

The remaining offset against the original closed-form number was
the closed form's own transport: at the edge state the gas is over
60% dissociated atoms, and CEA's mixture viscosity is 1.51x
Sutherland's pure-N2 value, which the form carries at the square
root (+23%) plus the Prandtl shift (+3%). With real-gas transport
the closed form reads 10.69 W/cm2 against the march's 10.7-11.2:
agreement to within a few percent. The 3-8% residual against a
frozen-Pr evaluation waits for the composition-aware wall closures.

Details and the honest caveats: [[numerics/foundations/fay-riddell|Fay-Riddell stagnation heating]].

## Dissociating nitrogen cylinder

The park5 chemistry gate (`cases/cylinder-m20-park5`): the Casseau
Mach-20 nitrogen cylinder, 1 m radius, pure N2 freestream at 6047 m/s
(M = 20.0), 400x200 grid, inviscid with slip walls, park 5-species
finite-rate chemistry, 3231 steps. The temperature and composition
fields and the stagnation-line profiles are on
[[numerics/foundations/five-species|the five-species page]].

| quantity | run value | published reference |
|---|---|---|
| shock standoff | 0.30 m | 0.25 m reacting, 0.29 m non-reacting (NS, 2-T) |
| stagnation-line peak T | 15.3 kK | 13.5-15 kK (2-T published); one-T ceiling 17.3 kK |
| stagnation-line N mole fraction at the wall | 0.34 | 0.1-0.25 (2-T published) |
| peak N mass fraction (near-wall layer) | 0.50 | order 0.1-0.25 mole |

The case names a Mach 20 nitrogen stream at 220 K, and the run flies that
state: the far field reads 220.0 K, 0.890 Pa, Mach 20.00. The gate in the
case directory holds it there and fails if the freestream drifts off the
declared state or the shock-layer scalars move more than a few percent.
Two faults used to bury that state. The uniform seed carried the air
mixture's energy into a pure-nitrogen march, and the supersonic inflow
handed the boundary pressure to the lane's internal-energy slot, which
pinned the far field to the 50 K floor and the whole family to a stream
nobody asked for. Both are fixed, and the shock-layer numbers moved less
than a kilokelvin, which says the enthalpy was always right and the
thermodynamics was not.

The dissociation is unmistakable: the frozen perfect-gas ceiling at
Mach 20 is 17.3 kK, and the reacting run peaks 2.0 kK below it with
the atomic-fraction field rising monotonically through the shock
layer. The one-temperature lane dissociates faster than Casseau's
two-temperature runs (their rates see a geometric-mean reaction
temperature), which is why my composition reads above their published
band while my stagnation temperature sits below their ceiling: the
chemistry converts more of the thermal pool into dissociation energy.

The standoff stays at the non-reacting value rather than shrinking to
the published 0.25 m. Two contributors I can quantify: this is an
inviscid run at 80k cells where the published numbers are thin-layer
Navier-Stokes at 155k cells with 2 um wall spacing, and the shock
rises over 2 cells so the crossing reads half a cell high. The
published standoff reduction is a shock-layer thickening effect
(recombination at the cool wall pulls the layer in); a slip wall
cannot reproduce it. The scalars that do transfer cleanly at this
resolution are the peak temperature and the composition field, and
both land in or above their expected bands for the right reason.

Drag from the pressure integral is not comparable at this resolution:
at the shoulder the wrapped shock sits on the surface cells, so the
integral picks up post-shock pressures and reads several times the
Newtonian 1.33. Surface-integrated drag from this lane is
measurement-limited, not citable. A body-fitted wall treatment is the
open step: the axisymmetric shell case already carries a cut-cell
implementation, and it tears the field down rather than closing the
gap, so the next move is a wall-local flux redistribution that keeps
the shock layer out of the average.

## Nozzle cold flow

Liou's CDV nozzle (AIAA 87-0355) as the axisymmetric Euler case,
at the archive's three exit pressure ratios. Deep-converged
marches (t = 160 to 400; the supersonic exit Mach is identical
between t = 160 and t = 400).

| case | regime | exit p/pt vs 1d | exit Mach vs 1d | mass flux spread |
|------|--------|-----------------|-----------------|------------------|
| pe/pt = 0.89 | subsonic | 0.874 vs 0.890 | 0.372 vs 0.411 | 15% plane mass flux spread (still settling) |
| pe/pt = 0.16 | fully supersonic | 0.197 vs 0.160 | 1.742 mass-avg vs 1.854 | 0.07% |
| pe/pt = 0.75 | normal shock | 0.729 vs 0.750 | shock at x ~ 6 vs 1d 7.56 | post-shock noisy |

The supersonic case is the cleanest gate: 5.8% centerline L2 against
the 1D profile, 0.07% mass flux spread across three planes, and an
exit Mach 6% under the 1D value because the staircase wall pinches
the effective throat area ~4% under the analytic contour. The subsonic
case holds its back pressure to 1.8% but trails the 1D exit Mach by
9.5% and still carries transients into its window. The shock case
captures the shock but stands it early with staircase noise downstream:
captured, not resolved, and the stair-stepped wall on the nozzle
remains the piece cut cells are meant to answer.

The continuum probe ships with these runs and every
other 2D/axisymmetric run: Bird's gradient-length local Knudsen
number per cell against the 0.05 breakdown threshold, in the run
log.

Details, the reservoir-inflow and back-pressure boundary
conditions, and the figures:
[[numerics/2d/nozzle-cold-flow|Nozzle cold flow]].

## The trim anchor

The trimmed capsule: 384x192, viscous, a 1600 K radiation-equilibrium wall,
1804 steps to t = 2 ms in 271.6 s. Peak surface pressure 0.959 MPa at
(x=1.549, y=1.281), maximum field temperature 25670 K, peak wall heat flux
0.575 MW/m2, C_D 0.501 with L/D 0.335. The peak pressure stays inside the
retired table's 1 MPa ceiling, so the march never needed the rebuilt
equilibrium-air grid to finish. This is the case the wake visualizations come
from, including the header render on the repository
front page.

## MPI scaling

The multi-rank Fortran march measured on the 34.6 km station case, 320x176,
30 steps, on my 8-core desktop (scaling plot on
[[numerics/parallel/mpi|the MPI page]]):

| ranks  | wall (s) | speedup |
|--------|----------|---------|
| serial | 1.5      | 1.00    |
| np=1   | 1.5      | 1.00    |
| np=2   | 0.8      | 1.88    |
| np=4   | 0.8      | 1.88    |

np=4 does not beat np=2 on this grid. At four ranks each column slab is 80
cells wide and the halo exchange plus the per-step body gather/broadcast cost
more than the thinner slabs save. The scaling page discusses the mechanism;
the honest summary is that the march scales to two ranks on entry-capsule
grids and needs wider slabs (a bigger grid, or a 3D decomposition) before
more ranks pay.

The parity discipline behind these numbers: every rank-count configuration is
gated against the serial Rust march, exact zero for the body-free cases and
4e-15 worst relative difference for the capsule body case.

## Shared-memory threading (prior result)

The thread-sharded 2D step, measured in the 1.5.0 era on a 384x384 mesh,
40 steps: 1.86x at four physical cores, 2.04x at six threads, backsliding to
1.75x at eight (hyperthread oversubscription). The step is memory-bound, so
the plateau tracks memory bandwidth rather than core count. Full table on the
[[numerics/parallel/threading|threading page]].

## Catalytic wall accommodation

The park5 viscous lane's wall-heating series
(`cases/cylinder-m20-park5-visc`): the same Mach 20 nitrogen
cylinder with `physics.mu = 1.716e-5` and an isothermal 1000 K
no-slip wall, with the wall accommodation stepped from the
non-catalytic wall to the fully catalytic one. All five runs march to
t = 3 ms and finish clean.

| accommodation | peak wall heat flux | against the non-catalytic wall |
|---|---|---|
| 0, the non-catalytic wall | 3.746 MW/m2 | 1.000 |
| 0.25 | 3.820 MW/m2 | 1.020 |
| 0.5 | 3.895 MW/m2 | 1.040 |
| 0.75 | 4.230 MW/m2 | 1.129 |
| 1, the fully catalytic wall | 4.893 MW/m2 | 1.306 |

Every leg puts the peak on the same face at (0.963, 0.315). The two
ends reproduce the committed non-catalytic and fully catalytic run logs
exactly, and the peak grows monotonically with the accommodation; the
gate in the case directory checks all three.

The absolute numbers stay withheld, because the flux comes from the
march itself and the flow at the wall is out of the range where a
continuum solution means anything. The continuum probe fails over the
stretch of wall that carries the peak, with the gradient length
Knudsen number up to 1089 and about 10 percent of the cells above the
breakdown threshold.

The earlier explanation for the gap between the walk-through peak and
the march's peak is withdrawn. It rested on two claims, and both were
wrong. The march does not hold a constant conductance of 4.93 W/m.K
across the field; the species lane mixes the RP-1232 species
conductivities with the Wilke rule, so its conductance varies with
temperature and composition. And the walk-through number was never a
model comparison: 0.240 MW/m2 came from a reporter outside the march,
on the stream the case never declared, 50 K at Mach 42 with an air
seed carried into a nitrogen march. Comparing the two was meaningless,
and the arithmetic that backed the claim is not reproducible.

An earlier note blamed the peak on a boundary cell carrying a whole
cell's volume across a straddled wall; sampling the row on four rays
and reading the cells next to the peak shows a smooth region instead,
16137 K above it and 18791 K below, so that explanation is withdrawn
too.

The mechanism the series was meant to show is real and visible in the
composition fields: the fully-catalytic run's atomic nitrogen collapses
toward the wall (0.09 against 0.02 in the first fluid ring) while the
hot shock-layer core is unchanged, and the added term
$$q_{cat} = \sum_s h^0_s \, J_{s,wall}$$
carries the 33.6 MJ/kg formation enthalpy of recombining nitrogen into
the wall. The accommodation coefficient spreads the walls that sit
between those two ends, so a measured surface efficiency maps onto the
quoted wall heat flux.

The flux comes from the march itself: the viscous pass evaluates the
wall flux on the same one-sided quadratic it applies to the field, and
the run log reports that number directly. See
[[research/claim-register|the claim register]] for the current state of
these figures.

Details: [[numerics/foundations/catalytic-walls|the catalytic
walls page]].

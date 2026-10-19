---
title: Research results
---

The numbers behind the solver features, measured on the cases I ship. Every
figure here regenerates from the case run logs in the repo, so the numbers age
with the code instead of drifting from it.

## Aeroshell station sweep

Four trajectory stations at 22.5 degrees angle of attack, equilibrium air,
viscous with a 300 K wall, 320x176 grid, marched to t = 15 ms. This is the
sweep that anchors the hypersonic 2D path; the 34.6 km station is the peak-q
point of my re-entry trajectory. March steps and wall time per station are
plotted on [[numerics/foundations/equilibrium-air|the equilibrium-air page]].

| station | steps | wall (s) | max density ratio | peak p (Pa) | max T (K) |
|---------|-------|----------|-------------------|-------------|-----------|
| 30.0 km | 1387  | 110.8    | 27.5              | 6.54e5      | 5149      |
| 34.6 km | 1976  | 160.7    | 35.7              | 9.63e5      | 7776      |
| 40.0 km | 2474  | 202.8    | 54.1              | 9.89e5      | 10271     |
| 42.5 km | 2632  | 219.4    | 67.0              | 1.09e6      | 14031     |

The march length now tracks the shock layer rather than a table ceiling.
Under the old table the two highest stations clamped at its 1 MPa top;
with the rebuilt grid (101 pressures to 10.1 MPa, 30 temperatures to
30000 K, sensible internal energy) the same stations read peak pressures
up to 1.09 MPa and the shock-layer temperatures rise accordingly. The
42.5 km station pays the most: its thin, hot shock layer needs the
smallest time steps, and its density ratio of 67 is the sternest test of
the march's positivity.

Lift-to-drag across the sweep runs 0.55, 0.61, 0.65, 0.72 from sea level
to 42.5 km: the cleaner freestream at altitude rewards the same 22.5
degree shell with a better ratio.

## Stagnation heating gate

The Fay-Riddell case (`cases/cylinder-stagnation-fay`): the Casseau
Mach-20 nitrogen freestream on a 1 m cylinder nose, 1000 K isothermal
wall, nose-only domain at global dt, 420k steps.

| quantity | value |
|---|---|
| edge plateau $\rho_e$, $T_e$, $p_e$ | 2.3e-4 kg/m3, 5040 K, 505 Pa |
| wall conductor, converged (refined 132x98 mesh) | 10.7 to 11.2 W/cm2 |
| conductor swing on the refined mesh | +-0.5 W/cm2 |
| closed form, pure-N2 Sutherland transport | 8.48 W/cm2 |
| closed form, CEA mixture transport at the edge state | 10.69 W/cm2 |
| agreement, converged march vs corrected closed form | 0 to 4% |

The 660k-step chain on the original 84x384 grid never settled: the
conductor held a 7.4-14.0 W/cm2 band indefinitely. The tangential
direction was the problem: the Hiemenz layer is ~11 mm and the
uniform 31 mm spacing aliased it to a third of a cell. The refined
132x98 grid (100 um wall cells, 2 mm at the stagnation line) kills
the breathing entirely: the conductor converges to 10.7-11.2 W/cm2
with +-0.5 swing and stays there across the settled window.

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
| stagnation-line peak T | 14.5 kK | 13.5-15 kK (2-T published); one-T ceiling 17.3 kK |
| stagnation-line N mole fraction at the wall | 0.38 | 0.1-0.25 (2-T published) |
| peak N mass fraction (near-wall layer) | 0.50 | order 0.1-0.25 mole |

The dissociation is unmistakable: the frozen perfect-gas ceiling at
Mach 20 is 17.3 kK, and the reacting run peaks 2.8 kK below it with
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
Newtonian 1.33. Cut cells are the honest fix for body-mounted
shocks, and until then surface-integrated drag from this lane is
measurement-limited, not citable.

## Nozzle cold flow

Liou's CDV nozzle (AIAA 87-0355) as the axisymmetric Euler case,
at the archive's three exit pressure ratios. Deep-converged
marches (t = 160 to 400; the supersonic exit Mach is identical
between t = 160 and t = 400).

| case | regime | exit p/pt vs 1d | exit Mach vs 1d | mass flux spread |
|------|--------|-----------------|-----------------|------------------|
| pe/pt = 0.89 | subsonic | 0.884 vs 0.890 | 0.414 vs 0.411 | +-5% (still settling) |
| pe/pt = 0.16 | fully supersonic | 0.192 vs 0.160 | 1.75 mass-avg vs 1.854 (peak 1.837) | 0.06% |
| pe/pt = 0.75 | normal shock | 0.755 vs 0.750 | shock at x ~ 6 vs 1d 7.56 | post-shock noisy |

The subsonic case is the cleanest gate (0.7% on both exit
quantities). The supersonic case chokes at the throat (M = 1.017)
and conserves mass flux to 0.06% across three planes; its exit
Mach sits 5.5% under the 1D value because the staircase wall
pinches the effective throat area ~4% under the analytic contour.
The shock case captures the shock but stands it early with
staircase noise downstream: captured, not resolved, the honest
weak point that cut cells exist to fix.

The continuum probe ships with these runs and every
other 2D/axisymmetric run: Bird's gradient-length local Knudsen
number per cell against the 0.05 breakdown threshold, in the run
log.

Details, the reservoir-inflow and back-pressure boundary
conditions, and the figures:
[[numerics/2d/nozzle-cold-flow|Nozzle cold flow]].

## The trim anchor

The finest shipped case is the trimmed capsule: 384x192, viscous, a 1600 K
radiation-equilibrium wall, 1772 steps to t = 2 ms in 239.1 s. Peak surface
pressure 2.07 MPa at (x=1.622, y=1.219), maximum field temperature 29379 K,
peak wall heat flux 12.2 MW/m2, C_D 0.501 with L/D 0.335. The peak pressure
sits far above the old table's ceiling, which used to clamp it at 1 MPa.
This is the case the wake visualizations come from, including the header
render on the repository front page.

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

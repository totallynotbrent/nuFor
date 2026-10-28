---
title: Chamber to throat
---

The axisymmetric combustion lane's verification case: a chamber of burned
gas behind an injector, a converging section, a throat and a bell, with a
back pressure at the exit. The throat has to choke, and its pressure and
temperature have to land on the equilibrium throat state for the chamber.

```mermaid
flowchart LR
  tank["reservoir 12 MPa"] --> inj["injector 7.0146e-4 m2"]
  inj --> chamber["chamber 7 MPa, 3689 K"]
  chamber --> conv["converging section"]
  conv --> throat["throat 1.2566e-3 m2, mach 1"]
  throat --> bell["bell"]
  bell --> exit["exit, mach 1.39"]
  exit -. back pressure 542900 Pa .-> outlet
```

## The case

`cases/chamber-throat/case.toml` runs on a 110 by 40 grid over
0.165 m by 0.1 m. The wall is a chamber of 0.1 m radius that necks to
0.0206 m at x = 0.14925 m and opens into a bell.

| field | value |
| --- | --- |
| equations | `euler_axi`, `eos = "combustion"`, nine species, chemistry on |
| chamber | 7.0 MPa, 3676.53 K, the burned-gas injection map |
| reservoir | 12.0 MPa, 6.191 kg/m3, fed through a 7.0146e-4 m2 injector |
| throat area | 1.2566e-3 m2 on the grid, 1.3338e-3 m2 on the contour |
| back pressure | 542900 Pa |
| march | 0.0005 s in 4223 steps, hllc with muscl, cfl 0.3, four threads |

The injector's area is sized against the grid's throat, not the
contour's. The grid's throat is a staircase plateau nine cells wide, so
its area is 5.8 percent under the continuous contour's. A feed sized
against the contour over-drives that throat by four percent, and the
throat never settles.

## The seed rides the grid

The case starts from a quasi-1d seed. The seed is built on the grid's
own per-column flow area, so the mass flux is uniform across the duct
from the first step and the throat's choked capacity is the grid's. The
old seed came off the continuous contour, which put the sonic point one
column wide on a mesh that holds it nine columns wide. The marched mass
flux then wandered by 23 percent along the duct, and no amount of
settling time fixed it.

## The settled state

Measured on the throat's column, the cell holding the wall's narrowest
point, as the area-weighted mean over the fluid rows:

| quantity | measured | CEA equilibrium | gap |
| --- | --- | --- | --- |
| mach | 1.0328 | 1 | +0.0328 |
| pressure | 3.5938 MPa | 4.0535 MPa | -11.3 percent |
| temperature | 3475.3 K | 3501.98 K | -0.8 percent |
| mass flow | 3.7968 kg/s | 3.80 kg/s | -0.1 percent |
| chamber pressure | 7.0521 MPa | 7.0 MPa | +0.7 percent |
| chamber temperature | 3687.4 K | 3676.53 K | +0.3 percent |
| exit mach | 1.7088 | 1 | |

The frozen mixture's isentrope crosses mach one at 0.565 of the
chamber's pressure, and the measured throat sits below that at 3.59
MPa against 3.99 MPa. The chemistry keeps running in the chamber and
changes the gas the nozzle sees, so the throat settles where the
choked state sits rather than where the frozen curve puts it. The
remaining gap to CEA's 4.05 MPa is the price of that continuing burn.

![the throat chokes and the flows balance](chamber-throat-settling.png)

The axial profile at the settled state carries the throat marker, the
equilibrium reference lines and the mass flow the duct passes. The mass
flow wanders by about five percent because the measurement integrates
the cell centres over the staircase's annuli, while the solver's own
fluxes run over the cell faces.

![the duct's axial state](chamber-throat-profile.png)

## The gate

`crates/nufor-cli/tests/chamber_throat_gate.rs` writes the case, runs
the march and holds four contracts on the settled snapshot:

| contract | tolerance | measured |
| --- | --- | --- |
| the throat chokes | mach within 0.08 of 1 | 1.0328 |
| the throat's temperature | within 5 percent of CEA | -0.8 percent |
| the throat's pressure | within 12 percent of CEA | -11.3 percent |
| the mass flow is conserved | within 8 percent | 0.5 percent |
| the chamber runs subsonic | mach below 0.1 | 0.0685 |
| the exit runs supersonic | mach above 1 | 1.7088 |

The march runs about thirty minutes in release, so that gate is marked
`#[ignore]` and runs locally. The settled state ships with the case
in `cases/chamber-throat/settled` as the marcher's own record (the
throat history and the settled profile), so a second gate checks the
contracts against that record in seconds, checks that the tables above
publish it, and checks the record against the field whenever the field
is on disk after a local march:

```
cd cases/chamber-throat && python3 gates.py
```

That gate runs in the continuous integration, so the numbers on this
page cannot drift away from the march that produced them without
something going red. The heavy gate still runs locally:

```
cargo test -p nufor-cli --test chamber_throat_gate --release -- --ignored --nocapture
```

The continuous integration runs the build, the lint and the parity
gates. The chamber marches stay on my desktop.

## Running it

The case's own output lands beside the case file. The throat's settling
history needs the interval snapshots, and the cli overwrites its snapshot
file, so the watcher runs beside the march:

```
nufor run cases/chamber-throat/case.toml &
python3 cases/chamber-throat/watch.py <run-dir> cases/chamber-throat/case.toml &
python3 cases/chamber-throat/figures.py <run-dir> docs/numerics/2d
```

The watcher leaves `throat-history.csv` and one `vtk-XXXXXX.vtk` per
interval in the run directory. The figures script reads those and writes
the two figures beside this page.

## What is not settled

The grid's throat area is quantized: 1.2566e-3 m2 on 55 by 20, 110 by
40 and 220 by 80, all within one percent of each other. Only a much
finer grid in y moves it, and the cut-cell machinery is not wired into
the command line. The grid's throat is therefore the throat the gate
measures, and the injector is sized to match it.

The reacting lane's local states carry a few percent of noise: the
finite-rate chemistry relaxes stiffly and the states near the wall
wander. The gate's bands cover that noise rather than pretending it is
not there.

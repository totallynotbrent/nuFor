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
  conv --> throat["throat 1.2566e-3 m2 on the grid, mach 1"]
  throat --> bell["bell"]
  bell --> exit["exit, flux mach 1.68"]
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
contour's. The grid's throat is a staircase plateau, so its area is 5.8
percent under the continuous contour's. A feed sized against the contour
over-drives that throat by four percent, and the throat never settles.

## The seed rides the grid

The case starts from a quasi-1d seed. The seed is built on the grid's
own per-column flow area, so the mass flux is uniform across the duct
from the first step and the throat's choked capacity is the grid's. The
old seed came off the continuous contour, which put the sonic point one
column wide on a mesh that holds it nine columns wide. The marched mass
flux then wandered by 23 percent along the duct, and no amount of
settling time fixed it.

## The settled state

The throat's pressure belongs to the choke, and the choke is found with
the flux mach: the mass-weighted mean mach, the same quantity the mass
flow balances on. It rises monotonically through the throat, so it gives
one crossing per run, and the crossing is interpolated between the two
columns that straddle mach one.

The mixed mach cannot do this. It squares the speeds and divides by the
reciprocals, and the staircase mixes the streams inside a cell, so at
the throat it swings wildly between neighbouring columns. Its crossing
wanders, and on the 220 by 80 grid it selected a column past the throat
altogether. Worse, it is not monotone in the way that matters: on the
shipped 110 by 40 case it crosses mach one three times, at 3.7131, 2.7677
and 2.5031 MPa. The figure this page used to publish was the first of
those three, read as if it were the only one.

Every number on this page comes from the flux mach. The values below are
the mean over the last eight samples of the march, with the spread those
eight samples carry:

| quantity | measured | band | CEA equilibrium | gap |
| --- | --- | --- | --- | --- |
| sonic pressure | 3.7483 MPa | 0.0258 MPa | 4.0535 MPa | -7.5 percent |
| sonic temperature | 3488.2903 K | 6.4 K | 3501.98 K | -0.4 percent |
| throat mass flow | 3.7892 kg/s | 0.0249 kg/s | 3.80 kg/s | -0.3 percent |
| chamber pressure | 7.0675 MPa | 0.0061 MPa | 7.0 MPa | +1.0 percent |
| chamber temperature | 3688.5862 K | 1.9 K | 3676.53 K | +0.3 percent |
| exit flux mach | 1.6769 | | 1 | |

The mass flow is the one number that lands on CEA almost exactly. The
temperature does too. The pressure does not, and the next section says
why.

![the throat chokes and the flows balance](chamber-throat-settling.png)

The axial profile at the settled state carries the sonic marker, the
equilibrium reference lines and the mass flow the duct passes. The mass
flow wanders by about five percent because the measurement integrates
the cell centres over the staircase's annuli, while the solver's own
fluxes run over the cell faces.

![the duct's axial state](chamber-throat-profile.png)

## What the gap is made of

A seven percent pressure gap can hide behind a hundred excuses, so the
audit takes it apart with the crate's own thermochemistry. The gas the
injector feeds is frozen by the time it reaches the duct: the
mass-flux-weighted composition holds across the whole converging
section, at 0.77 water, 0.11 hydroxyl and 0.07 oxygen, and the wall is
adiabatic. The gas that chokes is the gas the chamber state says it is,
and its properties come straight from the nasa-7 polynomials the crate
ships:

| property | value |
| --- | --- |
| mean molar mass | 18.559 g/mol |
| gas constant | 448.0 J/(kg K) |
| heat capacity at 3688.6 K | 3262 J/(kg K) |
| gamma | 1.159 |

Walking that gas's own frozen isentrope from the chamber's stagnation
state, with the mixture's own heat capacity integrated along the way,
puts the sonic pressure at 4.0422 MPa. CEA's equilibrium throat sits at
4.0535 MPa, which is 0.28 percent away from it. The gas is not the
explanation, and neither is the chemistry, and neither is the reference.

| reference | throat pressure |
| --- | --- |
| the marched gas's own isentrope | 4.0422 MPa |
| cea's equilibrium throat | 4.0535 MPa |
| the march | 3.7483 MPa |

What is left is the march itself. The converging section is audited on
its own from the last subsonic station, at x = 135.75 mm where the flux
runs at mach 0.61, and that station's own gas says the throat should
choke at 4.0955 MPa. The section is 8.48 percent short of its own
thermodynamics, and the shortage grows on the way in:

| x | marched pressure | off the isentrope |
| --- | --- | --- |
| 137.25 mm | 5.4404 MPa | -2.8 percent |
| 140.25 mm | 4.9461 MPa | -6.1 percent |
| 141.75 mm | 4.1905 MPa | -10.4 percent |
| 146.25 mm | 3.8121 MPa | -12.8 percent |

Nothing physical is taking the pressure out of that flow. The chemistry
has already finished, the wall does not conduct, and the gas's own
isentrope says the pressure should still be there. What is left is
numerical dissipation: the flow accelerates from mach 0.6 to mach 1
across thirteen columns and a sudden step in area, and the scheme bleeds
the transonic acceleration into entropy. Most of it happens where the
area steps into the plateau, at about 141 mm.

![the marched pressure against the gas's own isentrope](chamber-throat-isentrope-loss.png)

Knowing that the loss is ours changes what the fix is. The throat's
pressure is not seven percent low because the chemistry is late or the
closure is wrong. It is low because the mesh's staircased throat
dissipates the acceleration, and the cure is a mesh that represents the
throat's area properly, not a change to the gas.

## The wall's representation is the loss

The wall is drawn as a polygon and marched as a staircase: a cell is
fluid when its centre falls outside the body. That sounds harmless and
is not, because the throat is narrower than a grid row is tall.

On the 110 by 40 grid the polygon's neck sits at a radius of 20.605
millimetres and the rows are 2.5 millimetres thick, so the nearest row
boundary is at 20.000 millimetres. The staircase's throat is therefore
1.256637e-3 square metres against the contour's 1.333814e-3, a shortfall
of 5.79 percent.

The same row falls out on the 165 by 60 and the 220 by 80 grids, because
the contour's neck divides the same way against all three. The throat
does not improve with resolution. It stays wrong by the same amount.

Where the wall lands against the grid does matter, and that is the trap.
Scanning the row count found grids whose discrete throat area lands
within 0.13 percent of the contour by accident of alignment, and one
whose area is 2.2 percent over. Picking an accidental grid would have
posted a good number that depended on the wall's phase rather than on
the physics, so the fix has to make the throat's area independent of the
phase instead.

Holding the radial spacing fixed and doubling the axial resolution buys
almost nothing. Doubling the radial spacing instead buys two to four
points of deficit. That is the fingerprint of a radial representation
problem, so it settled what the fix has to be.

| grid | radial spacing | sonic pressure | the gas's isentrope | deficit |
| --- | --- | --- | --- | --- |
| 110 by 40 | 2.50 mm | 3.7483 MPa | 4.0422 MPa | -7.27 percent |
| 220 by 40 | 2.50 mm | 3.7275 MPa | 4.0341 MPa | -7.60 percent |
| 165 by 60 | 1.67 mm | 3.7972 MPa | 4.0364 MPa | -5.93 percent |
| 110 by 80 | 1.25 mm | 3.8189 MPa | 4.0258 MPa | -5.14 percent |

Two pairs of runs make the point. The 110 by 40 and the 220 by 40 share
a radial spacing and differ by 0.33 points of deficit. The 110 by 40 and
the 110 by 80 share an axial spacing and differ by 2.13 points. Axial
resolution is not what the throat needs; radial resolution is.

An earlier version of this table reported 3.7110, 3.7745 and 3.8397
megapascals, with deficits of 8.19, 6.49 and 4.39 percent, and it
claimed the deficit shrank on every level. Those numbers are withdrawn.
They came from the mixed mach, whose crossing wanders and, on the finest
grid, selected a column past the throat. With the choke located in the
flow the table above is the measured one, and it does not improve
monotonically.

Halving the Courant number does not change the picture either. At 110 by
40 with a CFL of 0.15 instead of 0.3 the sonic pressure lands near 3.68
megapascals, so the loss is not time-integration error.

Seeding the flow differently does not change it either. Seeding it from
the contour's true smooth areas instead of the mesh's quantised ones
changes the deficit by 0.1 percent. The march relaxes onto whatever
geometry it is discretised on, so the inception state washes out. The
loss is in the geometry, not in how the flow was started.

## The finest level is ill-conditioned

The 220 by 80 level is in the table the gate reads but it is not a
healthy run, and the gate excludes it from the trends rather than
quoting it as the best answer.

Its detected choke does not settle. Across its own snapshots the
pressure wanders between 3.40 and 3.56 megapascals and the crossing
alternates between two neighbouring columns:

| step | column | pressure |
| --- | --- | --- |
| 451 | 199 | 3.4154 MPa |
| 751 | 200 | 3.4063 MPa |
| 1251 | 199 | 3.4049 MPa |
| 1701 | 199 | 3.4431 MPa |
| 2451 | 199 | 3.5014 MPa |
| 3251 | 200 | 3.5373 MPa |
| 4701 | 200 | 3.4500 MPa |
| 5201 | 200 | 3.4681 MPa |

Under the old detector this level looked like the best result on the
case, at -4.39 percent, purely because the wandering crossing happened to
be read on a favourable column. It owed -13.16 percent once the choke
was located in the flow.

This is the staircase behaving as a staircase should. The finer the grid
gets near the throat, the smaller the individual cells are relative to
the corner they sit in, and the less the scheme can carry the
transonic acceleration across it. A representation whose throat area is
quantised to a whole radial row does not converge by being refined.

## The gate

`crates/nufor-cli/tests/chamber_throat_thermo.rs` reads the recorded
state, walks the gas's own isentrope with the crate's heats, and checks
the record against the recomputation. It writes the recomputed numbers
before it asserts on them, so a level that fails a contract still leaves
its measurements on disk to be diagnosed from.

```
cargo test -p nufor-cli --test chamber_throat_thermo
```

`cases/chamber-throat/gates.py` holds the contracts themselves and
checks that this page publishes the numbers the record carries:

```
cd cases/chamber-throat && python3 gates.py
```

| contract | tolerance | measured |
| --- | --- | --- |
| the throat chokes | the flux mach crosses 1 | it does, at x = 148.7 mm |
| the sonic temperature | within 5 percent of CEA | -0.4 percent |
| the sonic pressure | within 12 percent of CEA | -7.5 percent |
| the mass flow is conserved | within 8 percent | -0.3 percent |
| the chamber runs subsonic | mach below 0.1 | 0.002 |
| the exit runs supersonic | flux mach above 1 | 1.6769 |
| the gas agrees with CEA | within 1 percent | 0.28 percent |
| the march agrees with its gas | within 10 percent | -7.27 percent |
| the reference does not drift | the gas stays within 1 percent of CEA on every level | 0.93 percent |

The resolution study carries its own contracts, in
`cases/chamber-throat/resolution.py`. They describe what the
measurements have to show for the diagnosis to stand, rather than
wishing the grid were good: the gas has to keep agreeing with CEA, the
deficit has to be present and significant, the deficit has to track the
RADIAL spacing, and it must not track the axial one. A level that owes
more than the bound is ill-conditioned rather than informative, and it is
reported and excluded from the trends instead of being quoted as the
best level.

The shipped case and record are the 110 by 40 march, so the contracts
above describe that mesh. The grid summary in `settled/resolution.csv`
carries the finer grids.

Both gates run in the continuous integration. The heavy gate that
actually marches the case stays on my desktop, because it runs about
thirty minutes in release:

```
cargo test -p nufor-cli --test chamber_throat_gate --release -- --ignored --nocapture
```

## Running it

The march writes its snapshots beside the case file. The distiller walks
them, finds the choke in each and writes the record the gates and
the figures read:

```
nufor run cases/chamber-throat/case.toml
python3 cases/chamber-throat/distill.py <run-dir>
python3 cases/chamber-throat/figures.py cases/chamber-throat/settled docs/numerics/2d
```

The distiller writes `settled/history.csv` (the sampled sonic and
chamber states), `settled/profile.csv` (the axial profile),
`settled/stations.csv` and `settled/thermo_inputs.csv` (the gas along
the duct). The rust audit adds `settled/thermo_result.csv` and
`settled/isentrope.csv`, which is the curve the last figure draws. The
walked field itself stays out of the repository, and the gate checks the
record against it whenever a local march leaves it on disk.

A resolution study writes its record outside the case. The case and the
destination are explicit, so a finer grid cannot overwrite the shipped
110 by 40 record:

```
python3 cases/chamber-throat/distill.py <run-dir> \
  --case <run-dir>/case.toml --out <run-dir>/settled
```

`cases/chamber-throat/resolution.py` compares the distilled records
across grids and writes the summary the gate reads:

```
cargo test -p nufor-cli --test chamber_throat_thermo
python3 cases/chamber-throat/resolution.py --write \
  cases/chamber-throat/settled/resolution.csv <record dir> [...]
```

## What is not settled

The loss through the converging section is the open problem on this
case. Seven percent is the number the 110 by 40 mesh earns, and the
staircase does not refine its way out of it: the grid's throat area is
quantised at 1.256637e-3 square metres on the 55 by 20, 110 by 40, 165
by 60 and 220 by 80 grids alike.

The cure is to measure the throat's area from the body's own contour
rather than from the rows it lands on, so the answer stops depending on
where the wall falls against the grid. The repository already carries a
body-fitted axisymmetric grid for sphere-cone flows: its inner boundary
IS the body surface, it carries per-cell areas and per-face normals, and
it is mapped to the body instead of masked against it. A grid shaped like
that has a throat area that is right by construction and converges on the
contour's own figure instead of quantising to whole radial rows.

Adapting that machinery to this contour, to the nine-species chemistry
and to the viscous wall is the next piece of work, and the reference it
has to meet is the contour's own 1.333814e-3 square metres, independent
of the grid.

The reacting lane's local states carry a few percent of noise: the
finite-rate chemistry relaxes stiffly and the states near the wall
wander, so the sonic pressure moves by 0.7 percent over the settled
window. The gate's bands cover that noise rather than pretending it is
not there.

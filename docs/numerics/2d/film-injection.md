---
title: Blown-in film injection
---

An actively cooled wall carries two coolants, and they are not the same
physics. The regenerative coolant stays inside the wall's passages: it
takes energy from the wall and leaves by the chamber's own circuit, so
it changes the energy balance and not the chamber's mass balance. A
blown-in film crosses the wall into the gas. It brings mass, energy and
composition with it, so the chamber's mass balance has to carry it.

This page records the sealing of that crossing.

## What was wrong

The film's injection sat inside the wall's stated-temperature branch.
A wall that stated no temperature therefore injected nothing: the
coolant simply vanished, the chamber's mass balance forgot it, and the
leg the film was supposed to dominate looked like an adiabatic one.

The book's own counters did not help. They accumulated the film's area
and mass flow and were never read, so nothing checked the coolant's
amount and the two enumerations, one in the solver and one in the
report, had no reason to agree.

## What the seam guarantees now

The film's crossing is enumerated once, in `film_injection`, over the
wall faces the body's boundary crosses. It is independent of the wall's
thermal treatment: a coolant that crosses a wall has crossed a wall,
whether or not that wall radiates or convects.

The solver deposits that enumeration into the gas as mass, energy and
per-species mass. The wall report publishes the same enumeration as the
film's area and mass flow, so the number in the report is by
construction the mass the chamber took.

Where a boundary is cut obliquely, several faces can meet one gas cell.
The coolant's flux is then split across those faces by the same weights
the wall's own flux uses, so the film's mass lands where the wall's
flux lands.

## The gated result

`cases/film-injection/gates.py` builds and runs
`crates/nufor-core/examples/film_injection_gate.rs`, then writes
`cases/film-injection/gates.json`. The chamber is an eight-by-eight
annular grid with a flat wall at y = 0.02, the gas is resting
stoichiometric hydrogen and oxygen at seven megapascals, and the wall
states no temperature at all.

| quantity | measurement |
|---|---|
| film's mass flux | 5.000000 kg per square metre second |
| wall's stated temperature | 0 K, adiabatic |
| film's area | 0.010053 m² |
| reported coolant mass flow | 0.050265 kg per second |
| chamber's mass gain over one millisecond | 5.026548e-05 kg |
| mass the film books for that step | 5.026548e-05 kg |
| relative difference between the two | 1.213e-15 |
| nitrogen fraction in the receiving cell | 0.029762 |

The chamber's mass gain is positive under a wall with no thermal
treatment at all, it agrees with the reported coolant mass flow to
twelve decimal places, and the receiving cell's composition shows the
film's nitrogen.

The gate checks those four facts. The difference between the chamber's
gain and the film's booking is rounding noise fifteen orders of
magnitude below the coolant's mass.

## Reproduction

From the repository root:

```
python3 cases/film-injection/gates.py
```

The core behavioural tests are in `combustion_viscous_axi`:

```
cargo test -p nufor-core --lib combustion_viscous_axi
```

One test asks an adiabatic wall with a film to admit the coolant and
requires the chamber's mass gain to equal the film's booking. Another
requires the report and the solver to agree about the coolant's amount
and about the area it met.

Runs that sample the wall's duty now write two more columns to
`wall-duty.csv`: `film_area_m2` and `film_mass_flow_kg_s`.

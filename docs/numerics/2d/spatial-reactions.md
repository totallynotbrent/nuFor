---
title: Finite-rate reactions on a spatial composition
---

[[numerics/2d/spatial-chemistry|Spatial chemistry]] moves a scheme-width
composition through a flow. That is transport. Chemistry is the other
half: the reaction source that converts one parcel into another. This
page records coupling finite-rate GRI-Mech 3.0 kinetics to the spatial
composition, per cell.

## The seam

Each cell is a closed parcel for the duration of a source substep. The
chemistry a cell runs is exactly the zero-dimensional reactor's: the
same backward-Euler step on the fifty-three mass fractions, with the
temperature following from the cell's fixed internal energy. The
spatial pass and the parcel reactor share `backward_euler_step`, so
there is one solver, not two implementations that can drift.

`chem_react::apply_gri30_source` walks the field, hands each cell its
density and internal energy, and updates the composition in place. A
cell whose Newton solve does not converge is left at its incoming
composition and counted in the report, so a stiff cell pauses its
chemistry instead of corrupting the field.

No reaction subset is selected anywhere on this page. Every cell runs
the full mechanism. Adaptive selection is a separate question and
belongs to its own seam.

## The gated result

`cases/spatial-chemistry/reaction_gate.py` runs
`crates/nufor-core/examples/spatial_reaction_gate.rs` and writes
`cases/spatial-chemistry/reaction_gates.json`. Four methane-air parcels
from 900 K to 1700 K, lean to rich, are each advanced by one accepted
implicit step on both the spatial pass and the reactor's own step:

| quantity | measurement |
|---|---|
| cells reacted | 4 of 4 |
| cells deferred | 0 |
| worst Newton iterations | 8 |
| worst composition-sum error | 1.110e-16 |
| worst relative error vs the reactor | 7.153e-12 |
| worst absolute error vs the reactor | 5.183e-12 |

The field matches the parcel reactor to twelve significant digits. That
is rounding, not a discrepancy: both sides executed the identical step
function on the identical parcel.

## Reproduction

From the repository root:

```
python3 cases/spatial-chemistry/reaction_gate.py
```

The behavioural tests live in `chem_react`:

```
cargo test -p nufor-core --lib chem_react
```

One test requires a spatial cell to match the reactor's step on the same
parcel to 1e-9. Another requires a hot methane-air parcel to burn fuel
and form water while staying normalized. A third requires a pure
nitrogen cell to come through untouched, bit for bit.

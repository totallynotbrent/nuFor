---
title: Per-cell reaction-set selection
---

[[numerics/2d/spatial-reactions|Finite-rate reactions on a spatial
composition]] gave every cell the full 325-reaction GRI-Mech. Correct,
but wasteful: most species in a cell are trace passengers whose
evolution barely moves the parcel over one step. This page records the
selection seam that shrinks the evolved set per cell while keeping the
rate evaluation dense.

## What selection is, and is not

Selection chooses which SPECIES a cell's implicit step evolves. Frozen
species stay exactly at their incoming mass fractions for that step.
Crucially, freezing is a property of the state vector only: the rate
evaluation still sees all fifty-three species, so third bodies and
colliders keep their bath and reversible reactions keep their backward
channels. Only the solved degrees of freedom shrink. That separation
is what makes a reduced step comparable to the dense reference instead
of being a secretly different chemistry.

The mask is recomputed per cell from the incoming state on every pass,
ranked by `|dy/dt|` with a relative threshold and an absolute floor,
always protecting the major carriers (fuel, oxidiser, diluent, and the
stable products) and enforcing a minimum active count. A cell that
quenches sheds its radical set naturally on the next pass.

## The gated result

`cases/spatial-chemistry/selection_gate.py` runs
`crates/nufor-core/examples/spatial_selection_gate.rs` and writes
`cases/spatial-chemistry/selection_gates.json`. Five methane-air
parcels spanning 900 K to 1600 K and lean to moderately rich are
marched over 5 µs twice: once with the full mechanism, once with
selection. Both sides use the same adaptive substepping on refusal.

| quantity | measurement |
|---|---|
| mean evolved species per cell | 7.6 of 53 |
| evolved-set range | 7 to 8 species |
| steps refused (dense / selected) | 0 / 0 |
| worst major-species divergence vs dense | 7.063e-08 relative |
| worst temperature divergence vs dense | 7.260e-05 relative |

Divergence at this level is dominated by the frozen trace species'
feedback into the Jacobian conditioning, not by missing chemistry: the
majors track the dense trajectory essentially exactly.

## Honesty note

Seven active species is what the activity-ranked threshold earns at
these states. It is not a 53-species-accurate reduced mechanism in the
literature sense; it is a demonstrated bound on how many degrees of
freedom actually matter over a short horizon at these temperatures.
Tightening the threshold trades accuracy for cost predictably, and the
gate's tolerance fields are where that trade is asserted.

## Reproduction

From the repository root:

```
python3 cases/spatial-chemistry/selection_gate.py
```

Behavioural tests in `chem_select` prove a full mask reduces exactly to
the plain step, a partial mask freezes unselected species to machine
precision, and the protected majors always survive selection.

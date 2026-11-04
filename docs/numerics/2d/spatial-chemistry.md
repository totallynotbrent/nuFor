---
title: Scheme-width spatial chemistry
---

The combustion lane has two chemistry widths: a nine-species hydrogen
mechanism and GRI-Mech 3.0's 53-species methane mechanism. The generic
closure and the generic species transport already accept either width.
The missing link is the composition state that couples them inside a
spatial flow.

This page records the first sealed link in that chain: a conservation-
checked composition state whose species count follows its chemistry
scheme, advanced by the same conservative species transport used by the
nine-species lane.

## What the seam guarantees

`crates/nufor-core/src/chem_spatial.rs` provides a spatial composition
whose mass-fraction matrix has one row per scheme species. Construction
rejects a row whose length disagrees with the scheme, and rejects a row
that does not sum to one rather than normalizing it away. Normalization
after a transport step is permitted only for a positive, finite total
mass.

The advancement uses `species_transport::transport_species_generic`,
which takes the same mass-flux sweeps as the planar and annular flow
lanes and the same optional radial weighting for an annular grid. The
pass updates each species density and then reconstructs the mass
fraction against the new density. Every species row takes the same
route; none is translated through a nine-species adapter.

This is deliberately a transport and closure seam. It does not yet
couple GRI-Mech finite-rate chemistry to the spatial composition. In
particular, reducing the reaction set dynamically needs separate care:
a reversible reaction retains its partner, a broadened reaction retains
the species needed by its falloff form, and a third-body reaction retains
the collider contribution.

## The gated result

`cases/spatial-chemistry/gates.py` builds and runs
`crates/nufor-core/examples/spatial_chemistry_gate.rs`, then writes
`cases/spatial-chemistry/gates.json`. The measured field has 72 cells and
carries all 53 GRI-Mech species.

| quantity | measurement |
|---|---|
| methane mass fraction in the inlet column | 0.200000 |
| methane mass fraction initially in the receiving column | 0.005000 |
| methane mass fraction after one planar transport step | 0.200000 |
| transferred methane mass fraction | 0.195000 |
| least-positive species mass fraction anywhere in the field | 0.004500 |
| worst cell composition-sum error | 3.775e-15 |

The gate checks the scheme's 53-species count, that methane crossed the
step, that every species remained present, and that every cell stayed
within 1e-12 of a unit composition sum. The margin is more than twelve
orders of magnitude smaller than the composition scale.

## Reproduction

From the repository root:

```
python3 cases/spatial-chemistry/gates.py
```

The core behavioural tests are in `chem_spatial`:

```
cargo test -p nufor-core --lib chem_spatial
```

They cover a 53-species GRI-Mech composition, the generic thermodynamic
closure on that composition, a 65-row transport contact, and the original
nine-species scheme through the same interface.

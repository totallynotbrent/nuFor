---
title: case.toml
description: The case file — the reproducibility root of a nuFor run.
tags: [formats]
---
# case.toml

The case file is the reproducibility root of a nuFor run (spec 40). A case
directory holds `case.toml` plus, once a solver exists, `mesh/`, `results/`,
`logs/`, and `research/`. Every number that defines a run lives in this one
file, so a case can be archived, diffed, and re-run from scratch.

The schema is versioned: `schema_version` must equal the version the binary
understands (currently `1`). A case with a newer schema version is rejected
instead of being mis-parsed. Validation is implemented in
`crates/nufor-config` (Rust), which also exercises the example below as its
canonical test fixture.

## Full example

The tutorial case, `cases/tutorial/case.toml` (Sod shock tube):

```toml
schema_version = 1

[metadata]
name = "sod-shock-tube"
description = "Classic Sod shock tube on a coarse tutorial mesh"
case_revision = 1

[physics]
equations = "euler_1d"
gamma = 1.4
gas_constant = 287.0
unit_system = "si"
reference = { rho = 1.0, p = 1.0, u = 0.0, l = 1.0 }

[mesh]
nx = 100
x0 = 0.0
x1 = 1.0
source = "uniform"

[initial_condition]
type = "two_state"
left = { rho = 1.0, u = 0.0, p = 1.0 }
right = { rho = 0.125, u = 0.0, p = 0.1 }

[boundaries]
left = "wall"
right = "wall"

[numerics]
flux = "hll"
reconstruction = "first_order"
cfl = 0.5

[time]
final_time = 0.2
max_steps = 10000
residual_target = 1.0e-10

[output]
interval_steps = 100
formats = ["csv", "vtk"]
fields = ["rho", "u", "p"]
```

## Field reference

Unlisted keys are rejected: a misspelled setting must fail loudly, never
silently change a run.

| Key | Type | Required | Constraint | Meaning |
|-----|------|----------|------------|---------|
| `schema_version` | int | yes | = 1 | Schema version gate. |
| `metadata.name` | string | yes | non-empty | Case name, used in logs and results. |
| `metadata.description` | string | no | | Free text: goal, expected result. |
| `metadata.case_revision` | int | no (default 1) | >= 1 | Revision counter for the case definition. |
| `physics.equations` | string | yes | `euler_1d` | Equation set (2D sets later). |
| `physics.gamma` | float | yes | > 1 | Ratio of specific heats. |
| `physics.gas_constant` | float | yes | > 0 | Specific gas constant, J/(kg K) in SI. |
| `physics.unit_system` | string | no (default `si`) | `si` \| `nondim` | Unit convention for reported fields. |
| `physics.reference` | table | no | rho, p, l > 0 | Nondimensionalization reference state. |
| `mesh.nx` | int | yes | >= 2 | Number of cells. |
| `mesh.x0` | float | yes | | Left domain edge. |
| `mesh.x1` | float | yes | > x0 | Right domain edge. |
| `mesh.source` | string | no (default `uniform`) | | How the mesh is produced; a file path once external meshes exist. |
| `initial_condition.type` | string | yes | `uniform` \| `two_state` | IC shape. |
| `initial_condition` (uniform) | rho, u, p | yes | rho, p > 0 | Constant state. |
| `initial_condition` (two_state) | left, right | yes | rho, p > 0 per side | Left/right constant states, e.g. a shock tube. |
| `boundaries.left` | string | yes | `wall` \| `inflow` \| `outflow` \| `periodic` | Left boundary condition. |
| `boundaries.right` | string | yes | same set | Right boundary condition. |
| `numerics.flux` | string | yes | `hll` \| `hllc` | Riemann solver / approximate flux. |
| `numerics.reconstruction` | string | yes | `first_order` \| `muscl` | Spatial reconstruction. |
| `numerics.cfl` | float | yes | in (0, 1] | CFL number for the explicit time step. |
| `time.final_time` | float | no (default 0) | >= 0 | Stop time. |
| `time.max_steps` | int | no (default 0) | >= 0 | Hard cap on iterations. |
| `time.residual_target` | float | no | > 0 | Early stop when the residual drops below this. |
| `output.interval_steps` | int | yes | >= 1 | Write cadence in solver steps. |
| `output.formats` | string list | yes | non-empty, `csv` \| `vtk` \| `hdf5` | Writers to run. |
| `output.fields` | string list | yes | non-empty | Fields written, e.g. `rho`, `u`, `p`. |

## Validation rules

Parse errors (bad TOML, unknown keys, wrong types) reject the file
immediately. After a successful parse the following value rules are checked,
and every violation is reported at once:

- name non-empty; case_revision >= 1
- gamma > 1; gas_constant > 0; reference rho/p/l > 0 when present
- nx >= 2; x1 > x0
- every state has rho > 0 and p > 0
- CFL in (0, 1]
- final_time >= 0; residual_target > 0 when present
- at least one termination criterion: final_time > 0, max_steps > 0, or a
  residual_target
- interval_steps >= 1; formats and fields non-empty

## What is not stored here

Solver identity (binary version and code revision) is stamped at run time into
result and log provenance, not hand-edited into the case file. HDF5 holds the
large numerical datasets (mesh and results); `case.toml` holds the
interpretable description.

## Spec 40 mapping

| Spec 40 item | Where |
|------------|-------|
| Version the case schema | `schema_version` gate |
| Record solver version | run-time provenance stamp |
| Record code revision | run-time provenance stamp |
| Record numerical method selection | `numerics` section |
| Record physical units | `physics.unit_system` |
| Record reference conditions | `physics.reference` |
| Record mesh source | `mesh.source` |
| Record boundary-condition parameters | `boundaries` section |
| Record initial condition | `initial_condition` section |
| Record termination criteria | `time` section |
| Record output fields | `output` section |

## Related

- [[toml-config|Research — why `toml` parses this]]
- [[ffi-boundary|Research — the FFI boundary that consumes it]]
- [[architecture|Architecture]]
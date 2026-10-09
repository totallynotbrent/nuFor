---
title: immersed bodies and axisymmetric runs
---

Two case.toml additions extend the format beyond empty boxes and shock
tubes: an optional `[body]` section that puts an analytic solid in the
domain, and the `euler_axi` equation set for axisymmetric runs. This
page documents both. The base format lives in
[[case-toml]]; everything there still applies.

## The body section

`[body]` is optional. A case without it is an empty domain. With it,
the runner masks solid cells each step: any cell center inside the
body carries a quiescent reference state, and the flow solves around
it. The body is a body of revolution in axisymmetric runs and a
2D solid in planar runs.

The section is tagged by `type`:

```toml
[body]
type = "sphere_cone"
```

### sphere_cone

A spherical nose cap of radius `rn` tangent to a conical flank of
half-angle `delta_deg`, flaring out to base radius `rb`. This is the
classic capsule forebody shape. The sphere center sits at axial
station `xc`, so the nose tip lands at `xc - rn`; the tangent point
and base plane follow from the geometry.

| Key | Type | Constraint | Meaning |
|-----|------|------------|---------|
| `rn` | float | > 0 | Nose radius. |
| `delta_deg` | float | in (0, 90) | Cone half-angle, degrees. |
| `rb` | float | > 0 | Base radius. |
| `xc` | float | | Axial station of the sphere center; nose tip at `xc - rn`. |

### sphere

A circle: a cylinder in a planar 2D run, a sphere in an axisymmetric
run.

| Key | Type | Constraint | Meaning |
|-----|------|------------|---------|
| `r` | float | > 0 | Radius. |
| `cx` | float | | Center x. |
| `cy` | float | | Center y. |

### polygon

A closed polygon from its vertices, ordered around the boundary.

| Key | Type | Constraint | Meaning |
|-----|------|------------|---------|
| `verts` | [[x, y], ...] | >= 3 points | Vertex list in mesh units. |

## Axisymmetric runs

Set `physics.equations = "euler_axi"` for the axisymmetric Euler
equations: the planar fluxes with the annular update, meaning radial
faces weigh by their radii and the update carries the p/r geometric
source term.

Two conventions matter when writing such a case:

1. The mesh's y axis is the radius. `y0 = 0.0` puts the axis at the
   domain's south edge, which is the usual setup.
2. The bottom boundary is the symmetry axis, so it must be
   `slip_wall`: the mirror ghost makes the axis flux mass-free, which
   is the axis condition.

The wall mask from `[body]` composes with this: the capsule case
below runs the full immersed body at hypersonic freestream.

## Boundary kinds

The new kinds join the existing set (`wall`, `inflow`, `outflow`,
`periodic`, `profile_inflow`):

| Kind | Meaning |
|------|---------|
| `slip_wall` | Inviscid wall: the ghost mirrors the normal velocity. Also the axisymmetric symmetry axis. |
| `supersonic_inflow` | Fixed supersonic inflow state; see `[boundaries.inflow_state]` below. |
| `supersonic_outflow` | Supersonic outflow: the ghost copies the interior state, waves leave freely. |

A `supersonic_inflow` side reads its state from an adjacent
`[boundaries.inflow_state]` table rather than an inline value, so the
boundary list stays plain strings:

```toml
[boundaries]
left = "supersonic_inflow"

[boundaries.inflow_state]
rho = 1.0
u = 26.038
p = 1.0
```

`v` is optional and defaults to 0. Without an explicit
`inflow_state`, the inflow falls back to the uniform freestream from
the initial condition.

## The capsule case

`cases/capsule-peak-q/case.toml` ships as the worked example: a
sphere-cone capsule at a Mach 22 freestream, the peak-dynamic-pressure
trajectory point, run axisymmetric. The body's `xc = 3.9` places the
nose at x = 2.0 with the base plane near x = 7.1, inside the 9.5 m
domain.

```toml
schema_version = 1

[metadata]
name = "capsule-peak-q"
description = "sphere-cone capsule, mach 22 freestream, axisymmetric euler"
case_revision = 1

[physics]
equations = "euler_axi"
gamma = 1.4
gas_constant = 287.0

[mesh]
nx = 320
ny = 160
x0 = 0.0
x1 = 9.5
y0 = 0.0
y1 = 4.0
source = "uniform"

[initial_condition]
type = "uniform"
rho = 1.0
u = 26.038
p = 1.0

[boundaries]
left = "supersonic_inflow"
right = "supersonic_outflow"
bottom = "slip_wall"
top = "slip_wall"

[boundaries.inflow_state]
rho = 1.0
u = 26.038
p = 1.0

[numerics]
flux = "hllc"
reconstruction = "muscl"
cfl = 0.4

[time]
final_time = 0.5
max_steps = 200000

[output]
interval_steps = 5000
formats = ["vtk"]
fields = ["rho", "u", "p"]

[body]
type = "sphere_cone"
rn = 1.90
delta_deg = 9.5663
rb = 2.365
xc = 3.9
```

Run it with `nufor run cases/capsule-peak-q/case.toml`. The march
writes a VTK snapshot next to the case file: the bow shock forms
ahead of the cap with the post-shock layer hugging the surface, and
the wake trails the base.

## The browser workflow, upcoming

The case routes that make these files loadable, editable, and
runnable from the web UI (list, load, save, delete, run with streamed
frames) are the current work item. When they land, this page gains
the browser walk-through. The CLI path above works today.

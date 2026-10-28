---
title: Structured output
---

The 1D state has two plain-text writers, both driven from the Rust core
over the conserved fields:

- CSV: one row per cell (x, rho, m, e, u, p), ideal for spreadsheets and
  scripting.
- VTK: a legacy-format structured-points file (ascii) with rho, momentum, and
  energy as scalar arrays, readable by ParaView for snapshots.

The writers validate the array lengths up front and derive the primitives
(velocity, pressure) through the same kernels the solver uses.

## The 2D snapshot

A 2D run writes the same structured-points layout, one scalar field per
cell, and every field is available in every snapshot:

| Field | Meaning |
|-------|---------|
| `density` | the cell's density |
| `pressure` | the closure's pressure at the cell's state |
| `energy` | the conserved total energy |
| `temperature` | the closure's temperature |
| `mach` | the speed over the closure's sound speed |
| `velocity_x` | the axial velocity |
| `velocity_y` | the second velocity component: radial in the axisymmetric lane, transverse in the planar one |

A reacting run adds the composition, one field per species: `y_H2`,
`y_O2`, `y_H2O`, `y_OH`, `y_H`, `y_O`, `y_HO2`, `y_H2O2` and `y_N2` for
the combustion closure, `y_N2`, `y_O2`, `y_NO`, `y_N` and `y_O` for the
park five-species closure.

The velocity fields are what a reader needs for a mass flow. The mach
number carries the total speed, so a flow that turns has to be read from
the components to get the flux through a face.

## See also

- [[numerics/foundations/hdf5|HDF5 output]]
- [[numerics/foundations/restart|Restart]]
- [[numerics/foundations/diagnostics|Diagnostics]]
- [[formats/case-toml|case.toml format]]

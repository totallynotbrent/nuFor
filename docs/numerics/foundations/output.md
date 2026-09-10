---
title: Structured output
---

Two plain-text writers for the 1D state, both driven from the Rust core over the
conserved fields:

- CSV: one row per cell (x, rho, m, e, u, p), ideal for spreadsheets and
  scripting.
- VTK: a legacy-format structured-points file (ascii) with rho, momentum, and
  energy as scalar arrays, readable by ParaView for snapshots.

The writers validate the array lengths up front and derive the primitives
(velocity, pressure) through the same kernels the solver uses.

## See also

- [[numerics/foundations/hdf5||HDF5 output]]
- [[numerics/foundations/restart||Restart]]
- [[numerics/foundations/diagnostics||Diagnostics]]
- [[formats/case-toml|case.toml format]]

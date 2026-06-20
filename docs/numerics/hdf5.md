# HDF5 output

The solver can write its state to HDF5 through direct calls into the system
libhdf5 (no Rust binding), keeping the datasets the project owns:

- datasets: `centers`, `rho`, `m`, `e` (one per cell).
- the file reads back exactly what was written (round-trip tested).

The bundled script `tools/export_h5.py` reads a snapshot with h5py and prints a
CSV (x, rho, m, e, u, p) for plotting or further analysis on machines that have
Python and h5py installed.

## See also

- [[numerics/output|Output]]
- [[numerics/restart|Restart]]
- [[numerics/diagnostics|Diagnostics]]

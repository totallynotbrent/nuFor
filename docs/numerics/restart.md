# Restart

The run state can be saved to a compact binary restart file and loaded back
exactly: a header (magic bytes, cell count, gamma, time, step) followed by the
conserved arrays as little-endian f64. The load path validates the magic and the
block sizes, so a corrupt or foreign file is rejected rather than misread.

Round-trip tested: a saved state comes back bit-for-bit identical, including the
simulated time and the step count, which lets a run resume where it left off.

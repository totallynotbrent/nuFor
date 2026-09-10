---
title: Restart
---

The run state can be saved to a compact binary restart file and loaded back
exactly: a versioned header (magic, schema version, cell count, gamma, time,
step) followed by the conserved arrays as little-endian f64. The loader
validates the schema version, the exact file size, and the header physics, so
a corrupt file — or one written by a newer build — is rejected loudly rather
than silently misread.

Round-trip tested: a saved state comes back bit-for-bit identical, including the
simulated time and the step count, which lets a run resume where it left off.

## Compatibility

The on-disk layout carries an explicit schema version in the header. Current
files are written as version 1, and `nufor inspect` reports the version a file
was written with.

- **Legacy reads.** Files saved by the pre-versioned format (the earlier
  `NUFR1\0` magic) are still accepted and surfaced as version 0, so older
  result files keep working.
- **Forward rejection.** A file whose schema version is newer than this build
  understands fails with an explicit `IncompatibleVersion` error instead of
  being parsed at the wrong layout — old builds never guess at new files.
- **Strict shape.** The loader requires the exact expected byte count: a
  truncated file and one with trailing appended bytes are both rejected, so a
  half-written or tampered restart cannot come back as a plausible state.
- **Physical sanity.** gamma must be finite and greater than one and the time
  finite; a well-shaped file carrying NaN physics is refused.

This is the restart analogue of the schema policy: never corrupt a file to
force a read — reject, report, and migrate only when the change is deliberate
and documented.

## See also

- [[numerics/foundations/output||Output]]
- [[numerics/foundations/hdf5||HDF5 output]]
- [[numerics/foundations/diagnostics||Diagnostics]]
- [[formats/case-toml|case.toml format]]

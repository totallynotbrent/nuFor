---
title: Regression suite
---

The solvent-era rule that a solver stays trustworthy is simple: every time
anything changes, the whole thing must still produce the answers it produced
yesterday. nuFor's regression check is the automated test suite itself -- 29
named suites covering the 1D equation, the 2D Euler equation, viscosity, the
boundary conditions, the verification cases, the output writers, the CLI, and
the web server -- run as one gate.

## The gate

    tools/regression.sh

runs, in order:

1. `cargo fmt --check` -- formatting is part of the contract,
2. `cargo clippy --workspace --all-targets -- -D warnings` -- no warnings means
   no footguns,
3. `cargo test --workspace` -- every suite must pass,
4. a hygiene check -- a clean tree and no privately-tracked paths.

It exits non-zero on the first failing group, so CI and a developer both run the
same definition of "green". The GitHub Actions structure-check mirrors the
hygiene part; the native job mirrors the clippy + test part.

## What the suites anchor

The headline numbers the per-case tests lock down include:

- the Sod tube density field against the exact Riemann solution (L1 error
  shrinking as the mesh refines),
- the Lax shock tube staying positive and converging,
- the 2D HLLC cell matching the 1D reference along a plane,
- the oblique shock holding the exact post-shock state (the wedge),
- the two-shock regular reflection against its analytic triple,
- Poiseuille's parabola held by the viscous channel flow,
- the threaded 2D step being bit-identical to the serial one,
- the HDF5 and restart writers round-tripping exactly.

If any of those drifts, the regression gate turns red and the renumber is a
signal that a change moved a verified number, not that the test is fine to
patch around. (spec 189)

## See also

- [[numerics/verification|Verification ladder]]
- [[numerics/2d-euler|2D Euler solver]]
- [[numerics/3d-solver|3D solver]]

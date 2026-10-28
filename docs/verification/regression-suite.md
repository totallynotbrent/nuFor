---
title: Regression suite
---

The solvent-era rule that a solver stays trustworthy is simple: every time
anything changes, the whole thing must still produce the answers it produced
yesterday. nuFor's regression check is the automated test suite itself, run as
one gate.

## The gate

    tools/regression.sh

runs, in order:

1. `cargo fmt --check` -- formatting is part of the contract,
2. `cargo clippy --workspace --all-targets -- -D warnings` -- no warnings means
   no footguns,
3. `cargo test --workspace` -- every suite must pass,
4. the MPI parity gates, each in its own process so one gate's collective
   state cannot bleed into the next,
5. a hygiene check -- a clean tree and the path rules in
   `tools/repo_hygiene.py`.

Every group runs even when an earlier one fails, so one invocation reports
every problem in it rather than stopping at the first. The script exits
non-zero if any group failed.

The same rules run in CI. `.github/workflows/ci.yml` builds the workspace on
Linux and Windows, runs clippy and the test suite, and calls
`tools/repo_hygiene.py` for the path rules; the repo-hygiene job is a
seconds-long check with no CFD in it. A developer running the script and CI
running the workflow therefore agree on what "green" means.

Heavy gates stay local. Tests that march a real field for tens of thousands
of steps are marked `#[ignore]`, so they never consume CI minutes; run them
explicitly when the physics path they cover changes:

    cargo test --release -- --ignored

There are 17 such ignored gates spread across 11 test files: the combustor
ignition gate, the injector chamber gate, the chamber throat gate, the
equilibrium-air capsule artifact, the SA transport comparison, the viscous
validation lanes, and the convergence marches. Five of them are the MPI
parity gates, which the regression script runs one per process so their
collective state cannot bleed into each other.

## Counts

The suite is the workspace's own instrumentation, and the numbers move as
tests are added:

- 69 test binaries report success on a green run,
- 360 individual tests pass,
- 17 tests are ignored in the default run, and those are the 17 heavy gates
  above.

Reading those numbers off a run rather than quoting them keeps the page
honest as the suite grows.

## What the suites anchor

The headline numbers the per-case tests lock down include:

- the Sod tube density field against the exact Riemann solution (L1 error
  shrinking as the mesh refines),
- the Lax shock tube staying positive and converging,
- the 2D HLLC cell matching the 1D reference along a plane,
- the oblique shock holding the exact post-shock state (the wedge),
- the sphere-cone geometry reproducing the trade-study dimensions and
  the polygon body reproducing the cylinder bow shock,
- the axisymmetric sphere standoff at 1.2 percent of the
  Ambrosio-Wortman correlation, the case that proves the geometric
  source terms are present,
- the two-shock regular reflection against its analytic triple,
- Poiseuille's parabola held by the viscous channel flow,
- the 3D channel slab holding the same parabola across the span,
- the square duct holding its four-wall series solution,
- the 3D viscous operator draining a y-shear exactly like the z-shear,
- the threaded 2D step being bit-identical to the serial one,
- the HDF5 and restart writers round-tripping exactly.

If any of those drifts, the regression gate turns red and the renumber is a
signal that a change moved a verified number, not that the test is fine to
patch around.

## Case numbers

Tests cover the solver; the case campaigns have their own gates, and the
numbers they publish come from artifacts a fresh checkout can regenerate.
See [[research/results|Research results]] for the campaign tables and the
case each one comes from, and [[numerics/2d/nozzle-cold-flow|nozzle cold
flow]] for a worked example of a shipped config, the gate that reads its
output, and the tolerance the gate holds it to.

Three campaigns ship a gate beside their case today: the station sweep
(`cases/aeroshell-station/gates.py`), the dissociating nitrogen cylinder
(`cases/cylinder-m20-park5/gates.py`) and the stagnation chain
(`cases/cylinder-stagnation-fay/gates.py`). Each one reads the record of
its own run, measures the quantities its pages publish, and holds them to a
written tolerance. The stagnation gate reads the committed chunk table and
series behind it, and the dissociation gate reads `samples.csv`, the slices
`samples.py` distils out of the marched field, so both run in CI and a
drifting campaign number fails the build instead of ageing quietly in the
docs. The station gate reads the marched field itself, so it runs where the
campaign has been marched; the regression script skips a case that has not
been marched here and fails a case that has marched and drifted.

## See also

- [[verification/index|Verification ladder]]
- [[numerics/2d/euler|2D Euler solver]]
- [[numerics/3d/solver|3D solver]]

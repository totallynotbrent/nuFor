# Verification ladder

The scheme is checked rung by rung, each a stronger physical statement than the
last, so that a passing test means the solver is doing the right thing rather
than merely staying stable:

- Constant state: a flat field at rest stays exactly flat.
- Uniform advection: a uniformly moving field advects without distortion.
- Stationary shock: a normal shock placed in its own rest frame keeps its two
  plateaus and its jump.
- Sod and Lax shock tubes: compared against the exact 1D Riemann solution,
  with error norms that shrink as the mesh refines.

Each rung ships as a focused test with evidence; a successful run or a
screenshot is not verification.

## See also

- [[numerics/sod|Sod verification]]
- [[numerics/lax|Lax verification]]
- [[verification/gate-1d-euler|1D Euler gate]]
- [[verification/regression-suite|Regression suite]]

---
title: Diagnostics and log hardening
---

The solver reports *why* it stopped, and can scan a state for non-physical
cells, so a long run that blows up is diagnosed instead of silently producing
NaN.

## Termination reasons

Every run ends with one of these, exposed as `EulerResult.reason`:

- **Converged** — the residual dropped to or below `tol`.
- **TimeEnd** — the accumulated simulated time reached `t_end`.
- **MaxSteps** — the step budget (`max_steps`) ran out first.
- **BlowUp** — a non-physical cell appeared mid-run; the run stopped at that
  step so the result state is inspectable.

## Physical-validity scan

`check_physical(state, gamma)` scans every cell and reports whether all cells
have finite, positive density and derived positive pressure, plus the running
minima and the index of the first bad cell:

    PhysicalCheck { ok, min_rho, min_p, bad_cell }

The scan runs after every time step inside `euler_solve`, so a thermal or
mechanical overrun stops the run cleanly with reason `BlowUp` rather than
propagating garbage. Steady-state runs report `Converged` only when the residual
rule actually fires; a run exhausted by its step or time budget says so
explicitly rather than being misread as converged.

## See also

- [[numerics/positivity|Positivity]]
- [[numerics/output|Output]]
- [[numerics/restart|Restart]]

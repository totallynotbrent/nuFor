---
title: Lax shock tube
---

A second, much stronger Riemann problem used to stress the solver: larger
density and pressure jumps and a supersonic left state, producing a left
rarefaction, a central density discontinuity, and a strong right shock.

- Left: rho = 0.445, u = 0.698, p = 3.528. Right: rho = 0.5, u = 0.0, p = 0.571.
- gamma = 1.4, domain [-0.5, 0.5], time 0.15, transmissive boundaries.

The same wave family as Sod but with a supersonic left state and a much
stronger right shock, so the limiter and the positivity floor both get
exercised:

```mermaid
flowchart LR
    L["supersonic<br/>left state"] --> R1["rarefaction<br/>(stronger)"]
    R1 --> S["star region"]
    S --> C["contact<br/>(density spike)"]
    C --> Sh["strong right<br/>shock"]
    Sh --> Rt["right state"]
```

## Verification

The density is compared against the exact 1D Riemann solution, and the gate
holds three things (`crates/nufor-core/tests/lax.rs`): the density stays
positive through the strong wave, the momentum stays finite, and the density's
L1 error against the exact solution stays under 0.06 on a 200-cell mesh. A
second test marches 200 and 400 cells and requires the error to shrink, the
same first-order convergence the Sod case shows. Together with the Sod case
this forms the core of the shock-tube verification ladder.

## See also

- [[numerics/1d/sod|Sod verification]]
- [[verification/index|Verification ladder]]
- [[numerics/1d/euler|1D Euler]]

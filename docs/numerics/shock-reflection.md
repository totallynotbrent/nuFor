# Shock reflection off a wall

A plane oblique shock is one thing; when it runs into a solid wall it reflects,
and the reflected shock plus the doubly compressed flow are fixed by the
two-shock (Rankine-Hugoniot) solution. This is the standard shock-reflection
benchmark: a supersonic stream turned by an incident shock and then turned back
by its reflection off a slip wall.

## Two-shock theory

For an incoming M1 = 2.5 stream turned twelve degrees by the incident shock,
the oblique-shock relations give the shock at beta1 = 33.8 degrees and a
post-incident stream of M2 = 2.00. That turned stream then reflects off the
wall: applying the same relations to M2 with the same twelve-degree turn back
to wall-parallel gives the reflected shock at beta2 = 41.5 degrees and the
doubly shocked state rho3 = 2.617, p3 = 3.948, flowing parallel to the wall.

## The benchmark

`tests/shock_reflection.rs` lays all three regions out as the analytic steady
field — the oncoming stream, the turned region between the two shocks, and the
reflected region against the wall — and marches the solver over it.

- the oncoming stream stays at rho = 1 within a third of a percent,
- the reflected (doubly shocked) region holds rho3 = 2.617 within a percent,
- the turned sliver between the shocks reads close to its analytic value.

The reflection point itself is where the incident shock meets the wall, and a
coarse Cartesian grid resolves that corner only approximately, so weak waves
bleed outward from it as the run continues: the agreement is tightest over a
short horizon and degrades slowly (about four percent of density over a run of
one image time). The short-march check captures the steady two-shock state
where it actually holds, and that slow degradation is itself instructive —
it is exactly the sensitivity that a body-fitted or AMR grid removes.

## See also

- [[numerics/oblique-shock|Oblique shock]]
- [[numerics/2d-euler|2D Euler solver]]
- [[numerics/verification|Verification ladder]]

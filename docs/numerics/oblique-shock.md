# Supersonic wedge and oblique shock

A supersonic stream turning through a wedge is a canonical 2D compressible
problem. The sharper the wedge, the stronger the oblique shock that detaches
from it, and the whole shock is described exactly by the two-dimensional
theta-beta-M relations — a closed-form reference the solver is checked
against.

## The theta-beta-M relations

For an incoming mach number M1 and a turn through the wedge angle theta, the
shock lies at an angle beta from the incoming direction, found (weak branch)
from

    tan(theta) = 2 cot(beta) (M1^2 sin^2(beta) - 1) / (M1^2 (gamma + cos(2 beta)) + 2)

Once beta is known the downstream state is the normal-shock jump across the
normal component M1 sin(beta): density, pressure, and the downstream mach all
follow from the one-dimensional Rankine-Hugoniot relations.

## Verification

`tests/oblique_shock.rs` takes M1 = 2 and a 10-degree turn, computes beta and
the post-shock state from the relations, initializes a 2D domain split by a
shock line at that angle — the turned stream below it, the free stream above —
and marches with the solver. Because an oblique shock is a steady solution, the
solver must keep it in place:

- the free stream stays at rho = 1 exactly,
- the turned region holds the analytic post-shock density to about a quarter
  of a percent,
- the shock line does not drift.

Try the same exercise for a stronger turn (larger theta) and the post-shock
state follows the relations until the shock detaches.
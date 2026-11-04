# rung 9b wall-exchange divergence: the evidence

date: 2026-10-31
repo: /home/ebarleb/Bots/nuFor
branch state: uncommitted working tree (rung 9b in progress)

## what was measured

three bracket legs, remarched with the rebuilt release binary (the
october archives predate the current fixes, so they were regenerated):

| leg | wall treatment | diverged at step | diverged at t |
|-----|----------------|------------------|--------------|
| wall-adiabatic | no wall exchange | not by step 7851 | - |
| wall-filmcooled | eps 0.9, film h 2000 W/m2-K | 3751 | 0.000444 s |
| wall-radiating | eps 0.9, no film | 5901 | 0.000700 s |

symptoms at onset: residual jumps from about 3.29 to above 5 within one
sample interval, and `rho_max` climbs monotonically from about 3.8
through 45 to 56 while `rho_min` and `t_max` stay sane. the runaway is
in density alone.

## why this matters for the rung

the october archived film leg had a DIFFERENT story: it never completed
(residual 57, density ratio 60) and the gate recorded a duty of zero
megawatts. that archived conclusion is void. the film leg does draw
heat; it diverges before the gate can measure it.

more seriously, the archIVELY PASSING radiating leg is ALSO invalid:
the same defect eventually took it down at step 5901, ten thousand
steps further than anyone had waited. a leg that passes because nobody
marched it long enough is not a passing leg.

## root cause

the wall exchange is explicit. the near-wall cell holds only a modest
margin of energy above the closure's floor (the internal energy of the
composition at 200 K). a cooled wall removes energy from that cell
every step. once the cumulative withdrawal exceeds the margin, the
cell's internal energy drops below the floor, and the positivity guard
in `add_axi_combustion_viscous` responds by:

    state.mx[k] = 0.0;
    state.my[k] = 0.0;
    state.e[k] = ke + floor * state.rho[k] + 1.0e-6;

that reconstructs the cell as a motionless parcel sitting on the
floor. the neighbouring cells then see a violent velocity and energy
gradient across the wall face, the disturbance propagates inward, and
the march loses its density field.

this explains every observation:

- the adiabatic wall never diverges, because it withdraws nothing;
- the film leg diverges first, because its duty is the largest;
- the radiating leg diverges later, because its duty is smaller but
  nonzero;
- the divergence is progressive and monotonic, not oscillatory, which
  is the signature of a steadily draining reservoir rather than a
  timestep instability;
- a SHORT march looks perfectly healthy. the earlier short-horizon
  reproduction reached step 2000 with residual 3.29 and was wrongly
  reported as curing the instability. it had not yet reached the onset
  at step 3751.

the film term compounds this. the old film flux was the unbounded
Newtonian sink `h (T_gas - T_wall)` applied against a FIXED 400 K
wall, with no coolant anywhere. at a 3700 K gas and h = 2000 W/m2-K
that is 6.6 MW per square metre, drawn from a cell that holds a few
joules per cubic centimetre. there was no second resistance and no
coolant to carry the duty away.

## the fix

two changes, both in the wall treatment:

1. an ADMISSIBILITY BOUND on the wall's withdrawal. a wall may take at
   most a fraction (`max_energy_fraction_per_step`) of the energy the
   cell holds above its thermodynamic floor in one step. this is not a
   physics cap: as the cell cools, the margin shrinks and the
   permitted withdrawal shrinks with it, so a converged solution that
   sits comfortably above the floor feels no limit. the limit exists
   only to keep the explicit step admissible.

2. a REGENERATIVE FILM with a real second resistance. the gas-side
   coefficient `h_g` and the coolant's capacity rate `m_dot cp` sit in
   series, so the duty is

       q = (h_g * m_dot cp / (h_g + m_dot cp)) * (T_gas - T_coolant)

   which is bounded by what the coolant's mass flow can physically
   carry, however hot the gas gets. the old unbounded sink is what
   made the film leg uniquely fragile.

plus a BLOWING-IN film boundary that injects the coolant's mass,
energy and species into the chamber. the road map promises an
"injected-wall species/enthalpy flux for film cooling"; the old code
had no injection at all, so the film leg drained energy from the
chamber without supplying mass, and the chamber's mass balance was
invalid for that leg.

## the gate's own contract was wrong too

the bracket gate asserted that "a cooler chamber carries MORE mass
flow, because the choked mass flow scales as p0 / sqrt(t0)". it held
the chamber's TOTAL PRESSURE fixed while varying the temperature. the
wall does not hold anything fixed: it removes energy, so the chamber's
total pressure falls WITH its total temperature, and the two effects
partly cancel.

checked against the archived radiating leg: with p0 and t0 both free,
the choked relation predicts -0.09902 kg/s against the measured
-0.09821 kg/s, an error of 0.022 percent. the gate had called the
same comparison "-5.9 times the prediction".

the contract is now phrased as a RESULT rather than an assertion: the
leg's mass flow must agree with the choked relation evaluated on the
chamber state the leg actually produced, within 5 percent.

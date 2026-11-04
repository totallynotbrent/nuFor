---
title: GRI-Mech 3.0 chemistry
---

The methane/oxygen mechanism lane: the 53 species and 325 reactions of
GRI-Mech 3.0, with four-parameter Troe falloff, equilibrium-constant
reverse rates, and Cantera's third-body convention followed exactly.
This page describes what the mechanism is, how it was verified, and
where its tables live.

## What it is

GRI-Mech 3.0 is the standard methane/air mechanism: 53 species, 325
reactions, chemical-kinetics data drawn from the 2000 GRI consortium
compilation. Unlike the nine-species hydrogen lane, it carries carbon
chemistry, nitrogen chemistry, and heavy-hydrocarbon intermediates,
so it covers premixed flames, nonpremixed flames, extinction, and the
post-flame pollutant chemistry that a combustor's exhaust composition
depends on.

The mechanism is generated into Rust tables rather than transcribed by
hand. The generator walks the mechanism's own yaml with Cantera 3.2.0
installed, so every rate coefficient, efficiency table, and Troe
parameter comes out of a single audited source rather than out of a
retype. The generated file is checked in, so the build needs no Python
and no Cantera installation.

```
crates/nufor-core/tools/generate_gri30_tables.py    the generator
crates/nufor-core/src/gri30_generated.rs             the generated tables
crates/nufor-core/src/gri30_mechanism.rs             the mechanism's interface
crates/nufor-core/src/gri30_reversible.rs            rates of progress
crates/nufor-core/src/gri30_thermo.rs                NASA-7 polynomials
crates/nufor-core/src/troe.rs                       the falloff broadening
```

The generator pins the mechanism file's checksum before it reads
anything, so a silently swapped mechanism is an error rather than a
quiet change in every number downstream.

## The mechanism's own structure

The reaction set divides into three kinds, and the code keeps them
apart:

- plain Arrhenius reactions, with the modified form
  $k = A\,T^n\,\exp(-E/RT)$
- three-body reactions, whose rate picks up a collider bath
- pressure-dependent falloff reactions, whose rate blends between a
  high-pressure and a low-pressure limit through a broadening function

The falloff broadening follows Cantera's Troe form exactly. The central
factor carries two exponentials weighted by the mixture parameter $a$,
plus an optional third term:

$$
F_{\mathrm{cent}}(T) = (1-a)e^{-T/T_3} + a\,e^{-T/T_1}
  + \;\bigl[e^{-T_2/T}\ \text{if } T_2 \neq 0\bigr]
$$

The third term is added, not multiplied. The centre is then widened by
the Blanchard-Roscoe divisor, which is what puts the sharp bend in the
wings of the crossover where a lean flame's kinetics live. Writing
$w_0 = \log_{10} F_{\mathrm{cent}}$ and $l_{pr} = \log_{10} p_r$,

$$
c_c = -0.4 - 0.67 w_0, \qquad n_n = 0.75 - 1.27 w_0
$$

From those two coefficients the reduced pressure follows, and the
broadening factor with it,

$$
f_1 = \frac{l_{pr}+c_c}{n_n - 0.14(l_{pr}+c_c)}, \qquad F = 10^{\;w_0/(1+f_1^2)}
$$

and the rate constant itself is the blending

$$
\begin{aligned}
k &= k_{\infty}\frac{p_r}{1+p_r}F, &
\qquad p_r &= \frac{[M]\,k_0}{k_{\infty}}.
\end{aligned}
$$

A shorter parameter list is accepted by dropping only the terms it does
not use, and a reaction with no broadening falls back to the plain
Lindemann blend. Cantera's floor on both logarithms is honoured, so
the two implementations agree across the whole range rather than only
near the crossover.

The reversible reactions take their backward rate from the reaction's
equilibrium constant, which is built from the NASA-7 Gibbs free
functions of the participating species. Writing $\Delta\nu$ for the
reaction's change in molecule count and $p^{\circ}$ for the standard
pressure,

$$
K_p = \exp\!\left(-\frac{\Delta G^{\circ}(T)}{RT}\right),
\qquad
K_c = K_p\left(\frac{p^{\circ}}{RT}\right)^{\Delta\nu}
$$

and the natural logarithm of the concentration-form constant follows.
The Gibbs energies come from the same polynomial fits the enthalpy uses,
so the constant inherits the polynomials' accuracy instead of carrying a
second, independent fit. The backward rate constant is then derived from
the forward one through $K_c$, which keeps the two branches consistent
rather than transcribed.

The unit basis deserves care, because the constant and the densities it
pairs with have to agree. Cantera reports the concentration-form
constant on a mole-per-litre basis, and the mechanism's own rate laws
carry their concentrations in kilomoles per cubic metre. A kilometre of
moles spread over a cubic metre is the same amount of substance as a
kilogram-mole-scale number of moles in a litre, so the numerical value
of a concentration is the same in the two spellings and the constant
transfers without a factor.

For anyone who needs the concentration in moles per cubic metre instead,
the shift is $\Delta\nu \ln 1000$, and the implementation carries that
conversion as a separate function rather than folding it into the
constant. Doing it once, in one named place, keeps the two bases from
being confused.

Keeping the collider out of the stoichiometry matters here as well. The
collider cancels from both sides of the reaction, so counting it would
shift $\Delta\nu$ and bias every equilibrium constant that the reaction
touches.

## Third bodies and colliders

A three-body reaction proceeds at a rate proportional to the
concentration of a third body that carries away excess energy:
$$
q = [M]\left(k_f\prod_i [i]_{\text{excess}} \;-\; k_r\prod_j [j]_{\text{excess}}\right)
$$

There are two conventions for $[M]$, and mixing them up silently
changes the rate by a large factor, so this mechanism follows Cantera's
exactly.

For a generic `'M'` collider, every species can act as the third body.
Each species contributes its own concentration multiplied by an
efficiency, and a species the table does not mention contributes at the
default efficiency of one. The bath is that weighted sum:

$$
[M] = \sum_s \alpha_s\,[s]
$$

For a reaction with a *named* collider, spelled out in the equation
itself, only that named species can act as the third body. Its own
efficiency is one, and every other species contributes nothing:

$$
[M] = [s_{\mathrm{named}}]
$$

Concretely, `H + O2 + AR <=> HO2 + AR` uses argon alone as its
collider, whereas `H + O2 + M <=> HO2 + M` accepts every species with
its own efficiency. These are different reactions with different rates,
and they coexist in the same mechanism.

The equation above brackets both branches, so the same bath scales the
forward and the reverse term alike. That is not cosmetic. Were the
reverse branch to run without a collider, reversing the reaction would
move the collider instead of conserving it, and the detailed-balance
relation would break for every three-body reaction in the mechanism.
Equivalently, the backward rate constant is derived from the forward one
through the equilibrium constant, and the bath is then applied once to
the pair, so the same collider serves both directions.
The generated tables therefore carry, per reaction, a list of
third-body efficiencies and a default efficiency. The default is one
for a generic collider and zero for a named one, which makes the two
cases fall out of the same summation without a special case.

## Verification

The mechanism is checked against Cantera reaction by reaction, rather
than on a handful of summary numbers.

**Net rates of progress.** For every one of the 325 reactions, at five
temperature and pressure combinations ranging from 800 K and one
atmosphere to 3200 K and ten atmospheres, the rate of progress computed
from the Rust tables is compared against Cantera's own for the same
mixture. A fixture generated directly from Cantera holds the reference
values, and the test compares every reaction in the expanded list.

**Equilibrium constants.** The equilibrium constant computed from the
NASA-7 polynomials is compared against Cantera's for every sampled
reaction, which catches a transcription error in either the
thermodynamic coefficients or the stoichiometry.

**Falloff broadening.** The Troe broadening factors are compared
against Cantera's over a table of temperatures and reduced pressures
covering the mechanism's falloff reactions.


![the net rates of progress this mechanism computes against cantera's own](../img/gri30-net-rates-vs-cantera.png)

Each point is one reaction in one of the five evaluated mixtures. The
panel is drawn in the magnitude plane on logarithmic axes, because the
rates span forty orders of magnitude. Every point lies on the diagonal,
so the agreement is not confined to the fast reactions near the middle
of the range.

![the relative difference against cantera for every reaction and mixture](../img/gri30-cantera-agreement.png)

Most reactions agree to the last bit of a double, so their relative
difference is indistinguishable from zero and the axis floor stands in
for them. Twelve reactions rise above $10^{-14}$, and the worst of them
reaches a relative difference of $5.5\times10^{-13}$. Those are the
equilibria-sensitive reverse rates, where the smallest discrepancy in
the thermodynamic polynomials accumulates.

![the troe broadening factor against reduced pressure](../img/gri30-troe-broadening.png)

The mechanism mixes two kinds of falloff reaction. A plain one blends
its two limits smoothly, so its broadening factor is one by definition.
A Troe one carries four parameters, and its factor bends sharply away
from the crossover. The bend is what a lean flame's kinetics hinge on,
because it is exactly where the smooth approximation of the plain form
is furthest from the truth.

![the third-body bath for a generic collider and a named collider](../img/gri30-third-body-bath.png)

Both panels show the share of the bath that each species brings, taken
from the same mixture. Under a generic collider the bath is spread
across nine species according to their efficiencies, and the two majors
supply about three quarters of it between them. Under a named collider
one species supplies the whole bath, and every other species
contributes nothing at all.

**Third-body structure.** Separate tests pin the conventions
described above: a named collider draws only itself into the bath
while a generic collider lets every species join, and doubling the
bath doubles the net rate, which can only happen if both branches pick
up the same bath.

These tests are what turned up the largest defect in the mechanism's
bookkeeping. The third-body table had been reading the default
efficiency as one unconditionally, and the named-collider reactions had
been folding the collider into their own stoichiometry. Cantera's
convention differs from that assumption in exactly two places, and the
affected reactions include the hydrogen/oxidation branch that a
methane flame's earliest chemistry depends on.

## Marching a parcel

Rate-by-rate agreement is the mechanism checked against Cantera; burning
something is the whole kinetic chain, the thermodynamics, and the energy
coupling checked together. Both matter, because a mechanism can reproduce
every rate at a frozen state and still integrate to the wrong flame.

A zero-dimensional parcel makes that check possible. One kilogram of gas
in a rigid vessel, no losses, no geometry, so the only unknowns are the
chemistry and the heat release. Cantera has its own such reactor, so
there is an oracle for the whole trajectory rather than for one rate.

The case lives in `cases/chemistry/gri30-stoich-ignition/`:

```toml
[parcel]
mechanism = "gri30"
fuel = "CH4"
initial_temperature_k = 2000.0
pressure_pa = 101325.0
final_time_s = 2.0e-3
sample_interval_s = 2.0e-6
max_step_s = 5.0e-7
```

`nufor ignite cases/chemistry/gri30-stoich-ignition/case.toml` marches
it and writes the trace and the summary beside the case file. The initial
condition is stoichiometric methane and air at 2000 K and one atmosphere.

The march is not cheap per step. Each Newton iteration assembles the
Jacobian by forward differences, which costs one rate evaluation per
species per iteration, so a fifty-three species solve pays for fifty-three
evaluations each time. A run over a cool parcel's long induction period
therefore costs far more than the hot case quoted here, and the turn is
made by shortening the horizon, not by assuming the chemistry is slow.
What the solver buys for that cost is the agreement below: the forward and
reverse rates come from the mechanism as written, and nothing downstream
adjusts them to fit the answer.

At 2000 K the parcel's ignition delay is tens of microseconds, which
makes the whole burn fit inside a two-millisecond horizon. Cooler starts
behave the same way but take much longer: the same mixture does not
ignite at all at 1000 K even over a full second, and takes 43 ms at
1200 K. Cantera's own reactor gives these numbers, and
`tools/scan_ignition_delays.py` re-runs that scan against the mechanism
this repository pins, so the cooler-start claim is checkable rather than
remembered. Our own march agrees with Cantera there: the case in
`cases/chemistry/gri30-stoich-ignition-1200k/` ignites 43.251 ms after
the start, against Cantera's 43.35 ms, a difference of 0.23 %. The march
took 2,361,362 steps to cross that induction period, which is the cost
of the forward-difference Jacobian described above seen at its worst.
The shorter horizon is chosen for turnaround and not to flatter the
result.

![the parcel's temperature history against cantera's own reactor](../img/gri30-ignition-temperature.png)

The grey band is Cantera's trajectory and the orange line is the marched
one. Outside the flash the two coincide to within a few thousandths of a
kelvin. Through the flash the marched solution falls behind, reaching a
gap of 281 K before it catches up; the two panels below show what that
looks like and why it costs what it does.

![the flash window, enlarged](../img/gri30-ignition-flash-window.png)

The same two curves on their own axes over the first 0.12 ms. During the
flash the marched solution lags Cantera's by a noticeable margin, and
that lag is the price of a fixed backward Euler step through a stiff
event rather than an error in the mechanism.

![the major species histories against cantera's](../img/gri30-ignition-species.png)

Mass fractions of the five species that carry the burn, with Cantera's
drawn underneath in grey. Fuel and oxygen fall as carbon dioxide and
water climb, and the two sets of curves trace the same path.

![the temperature error through the burn](../img/gri30-ignition-error.png)

The error on a logarithmic axis, for the marched step and for one four
times coarser. Across the whole burn the median error of the marched
solution is 0.003 K, and only eighteen of its thousand sampled points
carry more than 10 K, all of them inside the first 64 microseconds.
Away from the flash the two solutions agree closely regardless of the
step, which says the settled state is a property of the chemistry rather
than of the integrator. Inside the flash the coarse run departs from
Cantera by more than a hundred kelvin while the fine run stays below the
same scale, so the error there is controlled by the step size and
shrinks when the step does.

The settled parcel is the strongest single number on this page:

| Quantity | nuFor | Cantera 3.2 | Difference |
| --- | --- | --- | --- |
| Final temperature | 3022.5617 K | 3022.5616 K | 0.0001 K |
| Final pressure | 165516.777 Pa | 165516.768 Pa | 0.010 Pa |

A tenth of a milliKelvin on the equilibrium temperature of a methane
flame, on a run with no fitted constants anywhere in the chain. The
comparison takes both reactors at the same instant, the end of the
marched trajectory, rather than comparing their own last samples.

**How the march works.** The parcel is advanced with backward Euler and a
Newton solve for the mass fractions at each step. Two details decide
whether such a march converges at all. The Jacobian is scaled by the step
size, because the residual is the state minus the state's evolution over
that step; leaving the step out of it makes the solve ill-conditioned
exactly where the stiffness is greatest. And the finite differences are
taken relatively, skipping species that carry no measurable mass, since
perturbing a vanishing radical by an absolute epsilon fabricates a
radical that the chemistry then reacts with. Neither is exotic, and both
were found by watching a convergent-looking march refuse to light.

## What is not wired yet

The mechanism's tables, thermodynamics, rates, and falloff broadening
are complete and verified against Cantera, and it now burns a parcel from
a case file with the settled state matching Cantera to a milliKelvin.

Two things are still ahead of a spatial case selecting it. The
nine-species hydrogen lane still owns the combustion closure, and it
carries a fixed species count through its state vectors, its closure, and
its transport, so widening that to fifty-three species touches roughly
170 call sites. And the marched parcel is a zero-dimensional reactor,
where no transport, no grid, and no boundary condition exist. Putting
methane into a flame means solving the mixture and its transport on a
mesh, and that is the harder half of the work.

The separable groundwork for both is in place. The transport routines
take slices and a species count rather than a fixed nine-species
layout, and they reproduce the legacy nine-species table exactly, six
tests over six properties comparing bit for bit. The thermodynamic
closure likewise works over a slice, and against the hydrogen lane it
matches bit for bit as well. The general scheme interface selects
between the two mechanisms without a branch at every call site.

So a case can select `gri30` and march a parcel today. A flame on a
grid still cannot.

## Where the numbers live

The mechanism tables and the falloff oracle come from
`crates/nufor-core/tools/generate_gri30_tables.py`, the thermodynamic
coefficients from `tools/gen_gri30_thermo.py`, the transport table from
`tools/gen_gri30_transport.py`, and the rates of progress live in
`crates/nufor-core/src/gri30_reversible.rs`. The fixtures holding
Cantera's reference values are generated by `tools/gen_gri30_netrates.py`
and `tools/gen_gri30_eqfixture.py`.

The reaction verification figures on this page are drawn by
`cases/chemistry/figures.py` from the numbers that
`crates/nufor-core/examples/gri30_verification.rs` prints, so they show
the crate's own evaluation rather than a re-derivation in a scripting
language.

The ignition figures are drawn the same way, by
`cases/chemistry/ignition_figures.py`, from the trace the crate wrote
when the case was marched and from the oracle trace that
`tools/gen_ignition_reference.py` produced with Cantera's own reactor.
Neither is transcribed by hand.

Both families rebuild together with one command:

```
bash cases/chemistry/regen-figures.sh
```

It regenerates the crate's verification csv, replots the mechanism
figures and the ignition figures, and writes them where the page reads
them. The pictures come out byte for byte the same as the ones published
here.

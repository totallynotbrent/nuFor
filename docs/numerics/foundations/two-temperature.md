---
title: Two-temperature relaxation (Park 2T)
---

The equilibrium-air closure assumes every energy mode shares one
temperature. Behind a strong shock that is wrong on the microsecond
scale: the translational pool jumps across the shock while the
vibrational modes fill over many collisions. The two-temperature rung
carries that lag as a live field.

## The model

One vibrational energy per cell, $e_v$ in J/kg, for the nitrogen
molecule (harmonic oscillator, characteristic temperature 3393 K). The
field is advected with the Euler mass fluxes (donor-cell) and relaxed
toward the translational pool with the Landau-Teller source

$$\frac{d e_v}{d t} = \frac{e_{eq}(T) - e_v}{\tau}$$

integrated backward-Euler each step,

$$e_v^{new} = \frac{e_v + (\Delta t/\tau)\, e_{eq}(T)}{1 + \Delta t/\tau}$$

$\tau$ is the Millikan-White correlation with Park's high-temperature
collision-limited cutoff, evaluated per cell from the local translational
temperature, pressure, and number density. Backward-Euler is
unconditionally stable and exact in both limits: frozen flow at small
$\Delta t$, equilibrium at large $\Delta t$.

The vibrational pool lives inside the conserved total energy: the
translational temperature reads `ein_tr = ein_full - e_v`, so the
relaxation moves energy between pools without touching the conserved
sum. Total energy is conserved by construction, which the gate checks.

## Closure and seeding

The Euler fluxes stay perfect-gas nitrogen (gamma and R from the
case). The vibrational field seeds at the initial condition's
equilibrium value, so a cold freestream starts vibrationally cold and
the shock layer is where the split develops. The reaction temperature
coupling (`sqrt(T * T_v)`, the Park mixing rule for dissociation
rates) is already implemented in the rates module and activates when
the species-continuity rung lands.

## Validation

Three gates pin the rung:

- the `T_v(e_v)` inversion round-trips the caloric form to 1e-9
- a thermostatted box at fixed translational temperature tracks the
  analytic Landau-Teller decay curve
  $e_v(t) = e_{eq}(T)\left[1 - e^{-t/\tau}\right]$ over three
  relaxation times
- a coupled box (temperature free to drift) conserves energy to
  machine precision and converges to the analytic coupled fixed point
  `e_v = e_vib(T(ein0 - e_v))` to 0.1%

The external anchor is the Mach 20 nitrogen cylinder of Casseau et
al. (2016): freestream 6047 m/s, 0.89 Pa, 220 K over a 1 m radius
cylinder. The inviscid lane reproduces the classical post-shock
translational temperature (17-18 kK on the stagnation line) with a
clear vibrational lag (T_v near 12-13 kK, a T/T_v split of ~1.4). The
quantitative standoff comparison against the published hy2Foam and
dsmcFoam profiles needs the wall-resolved viscous lane, which is the
Fay-Riddell gate's work (see the roadmap ladder).

![park 2T gate curves](../research/park2t-gates.png)

![casseau stagnation line](../research/casseau-stagnation.png)

## Using it

Set `physics.eos = "park2t"` on an `euler_2d` case:

```toml
[physics]
equations = "euler_2d"
gamma = 1.4
gas_constant = 296.8
eos = "park2t"
```

The run writes the vibrational temperature field alongside the main
snapshot as `<case>-tv.vtk` (scalar `t_vib`, K), so the T/T_v split
can be plotted directly against the translational temperature field.

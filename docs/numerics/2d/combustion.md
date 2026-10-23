---
title: H2/O2 combustion
---

The nine-species hydrogen/oxygen lane: a finite-rate mechanism coupled
to the Euler march through a species field, with a closure that turns
the conserved state into face pressures and sound speeds. This page
covers what the mechanism is, how the coupling works, and what its
numbers were checked against.

## What it is

The `combustion` closure models a reacting hydrogen/oxygen mixture with
the nine species H2, O2, H2O, OH, H, O, HO2, H2O2 and N2, evolved by
the 19-reaction mechanism of Li et al. (2004). Where the park5 lane
models air dissociation for reentry flows, this lane models a rocket
combustor: a fuel-oxidizer mix that enters unburned, ignites on its
own reaction heat, and relaxes toward equilibrium while it advects.

A case selects the lane with `eos = "combustion"` under
`equations = "euler_2d"`. The initial fill is stoichiometric unburned
H2/O2 by mass (2.016/18.016 H2, 16.0/18.016 O2), set by the march. The
reaction terms can be switched off with `chemistry = "off"` in the
`[physics]` section, which freezes the composition and turns the lane
into a multispecies transport test of the same machinery.

## The mechanism

The reaction set is the 19-reaction H2/O2 mechanism of Li, Zhao,
Kazakov and Dryer (2004): seven H2/O2 chain reactions, four
dissociation-recombination reactions and the HO2/H2O2 peninsulas. Each
forward rate is the modified Arrhenius form

$$k_f = A\,T^n \exp(-E / R T),$$

with the A, n, E triplets transcribed from the mechanism table in SI
units, third-body efficiencies folded in, and duplicate reactions
summed. The reverse rates come over each reaction's equilibrium
constant from the same NASA-7 polynomial fits that supply the
thermodynamics, so forward and reverse are consistent by construction.

Each species' enthalpy and specific heat are the NASA-7 polynomials,
the standard 200 K to 6000 K curve fits with the formation enthalpy
baked into the low-temperature intercept. The fits clamp below 200 K,
and the closure matches that clamp exactly: an energy inversion that
ignores the clamp turns non-monotone below the floor and the bisection
walks off its bracket.

## The closure

The conserved state carries total energy; the march subtracts the
kinetic part and hands the internal energy and the nine mass
fractions to the closure. The closure answers the three quantities
every face needs:

$$p = \rho R_{mix} T, \qquad a = \sqrt{\gamma_s p / \rho}, \qquad
\gamma_s = c_p / (c_p - R_{mix})$$

with $R_{mix}$ the mass-weighted mixture gas constant and $c_p$ the
mass-weighted NASA-7 specific heat, both evaluated at the temperature
solved from the energy. That temperature inversion is a Newton
iteration on the energy curve with the analytic derivative (the
mixture $c_v$), falling back to bisection inside the 200 K to 6000 K
fit window when Newton stalls. The frozen sound speed (composition
held fixed through the wave) is the honest acoustic closure for a
finite-rate gas: equilibrium sound speed would need the composition
derivative through the mechanism, which the frozen form deliberately
avoids.

The internal energy carries the formation offset, so an unburned
fill's specific internal energy at 298 K sits near minus 206 kJ/kg
relative to the elements' zero. That offset is exactly the reaction
heat: burning the mix at fixed energy converts the offset into
temperature, and nothing else in the march needs to know a flame
exists.

## The march

The step mirrors the park5 lane: a Heun two-stage update of the
conserved state with the species transported on the same face mass
fluxes, then the reaction source as an operator-split substep at the
post-step state. The species advection is the same upwind-on-HLLC
transport the park lane uses, generalized to any species count. The
chemistry substep adapts its internal substep to the fastest reaction
timescale and caps the substep count per step, which keeps the stiff
ignition transient inside the step budget.

The transport is honest about one limitation: the near-face cell feeds
itself back through an inflow unless the case says otherwise. The
`[physics.injection]` block names the species mass fractions the case
injects, and an inflow face carries those fractions instead of the
near-face cell. The block also sets the initial state, so a chamber
running on burned gas starts at that composition's own energy rather
than the premix's. A case with no injection block runs the
stoichiometric unburned fill.

## Verification

Three gates run in the library against NASA CEA 2024 reference points,
banked from real CEA runs of the same reactant basis (1 mol H2 + 0.5
mol O2 at 70 bar):

1. The flame point. Seeded with CEA's equilibrium composition at its
   adiabatic flame temperature (3676.53 K at 70 bar), a short reaction
   burst barely moves the composition and the closure round-trips the
   flame energy to the anchor temperature within 15 K. An explicit
   cold-start induction at 800 K cannot be resolved in test time, so
   the equilibrium-hold form is the honest gate.
2. The equilibrium composition at fixed (T, p). From the unburned
   fill at 70 bar, the mechanism relaxes to its fixed point and the
   mass fractions compare against CEA's at 3000 K and 3400 K. The
   water branch carries the known mechanism-versus-CEA gap (the
   formation-data vintage of the OH fit is the dominant term), so the
   tolerance there is 0.10 in mass fraction and the radical branches
   get a looser band.
3. The chamber contract. At CEA's chamber state (70 bar, 3676.53 K,
   3.612 kg/m3) the closure reproduces the banked chamber pressure
   within 2 percent and the sound speed within 5 percent of CEA's
   1480.49 m/s, which is what a throat expansion consumes downstream.

The case-level gate is the combustor duct: a planar constant-area
channel with the unburned mix entering at chamber conditions (70 bar,
800 K), reacting through its own heat release, and the steady exit
state compared against the CEA equilibrium at the same station.

## Walls and transport

With `physics.mu` set, the lane carries molecular transport: momentum
shear, heat conduction, and species diffusion. The per-species data is
the GRI-Mech 3.0 transport database (Lennard-Jones pairs), viscosity
by the Chapman-Enskog form with the Neufeld collision integrals,
conductivity by the Eucken correlation, and mixtures by the Wilke
rules, the same mixture-averaged fidelity class as the air lane's
RP-1232 fits. Binary diffusivities use the classic Bird-Stewart-
Lightfoot form; the H2 pair's known 35-50 percent overprediction of
light-gas diffusivities is a documented property of the database pair
and is carried, not patched.

A no-slip isothermal wall (`wall_temperature` set) carries one-sided
quadratic shear, the Fourier conduction at the wall temperature, and
zero-gradient species (a non-catalytic wall). Two limits are honest
here: the wall flux sits at the resolution of the first cell, so a
thermal boundary layer much thinner than the near-wall cell barely
moves the mean field (the cooled-wall bracket at 70 bar with 5 mm
cells cools the first cell by degrees per flow-through, not hundreds
of kelvin), and the film-cooling enthalpy-injection wall rides the
same species-plus-enthalpy flux machinery the transport already
carries, with no separate model.

## Where the numbers live

The mechanism and fits are `crates/nufor-core/src/combustion.rs` and
`combustion_thermo.rs`, the closure `combustion_closure.rs`, the march
`combustion_march.rs`, the species transport shared with park5 in
`species_transport.rs`, and the CEA gates in
`combustion_cea_gates.rs`.

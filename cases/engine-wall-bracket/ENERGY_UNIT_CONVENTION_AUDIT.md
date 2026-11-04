# Energy unit convention audit

## What was found

The axi engine lane's cooled-wall legs diverged: the film leg blew up at step
3751 and the radiating leg at step 5901, both with `rho_max` climbing from
about 3.8 to 45 and then 56 while the residual, `rho_min` and `t_max` stayed
sane. The adiabatic leg survived past step 7851.

The immediate suspicion was a numerically unstable wall exchange. That was
wrong. The actual cause is a unit error in how the solver reads its own
conserved state, and it sits in three places.

## The convention

`crates/nufor-core/src/state2d.rs` defines the conserved state:

```rust
/// total energy rho*e_t per cell.
pub e: Vec<f64>,
```

and `prim_to_cons2d` builds it as `e[i] = r * et[i]`, where `et` is the
SPECIFIC total energy in J per kilogram. So `state.e` is an energy DENSITY in
J per cubic metre, equal to `rho * e_total`.

Meanwhile the closure's floor, `combustion_closure::internal_energy_at`,
returns a SPECIFIC internal energy:

```rust
let h = ct::mix_h_mass(y, t);   // J per kilogram
h - mix_gas_constant(y) * t     // J per kilogram
```

Measured on the axi lane's own fixture (rho = 12.64 kg/m³):

```
state.e                       = 8.930158e+6   J/m³   (energy density)
state.e / rho                 = 7.064998e+5   J/kg   (specific total energy)
specific kinetic energy       = 6.320000e+4   J/kg
closure floor at 200 K        = -3.720868e+5  J/kg
```

## The three defects

### 1. the positivity guards were inert

Five guards across four modules subtracted a SPECIFIC kinetic energy from the
energy DENSITY and compared the result against a SPECIFIC floor:

```rust
let ke = 0.5 * (state.mx[k]^2 + state.my[k]^2) / rho;   // J/kg
let e_int = state.e[k] - ke;                            // J/m³ - J/kg
if ... || e_int < e_floor[k] { ... }                     // J/m³ < J/kg
```

The mixed quantity came out around 8.87e6, while the floor was -3.72e5, so
the comparison never tripped on a healthy cell. The guard was watching for
something that could not happen.

Its restore path was worse. It mixed the conventions in reverse:

```rust
state.e[k] = ke + e_floor[k] * state.rho[k] + 1e-6;   // J/kg + J/m³
```

and, in the axi lane, it also dropped the cell's momentum to zero, which is
what turned one cooked cell into a poisoned neighbour and a runaway.

Sites fixed: `combustion_march.rs` (two), `axi.rs` (three),
`cutcell_combustion.rs` (one), `solver_u.rs` (two),
`combustion_viscous_axi.rs` (one). Every one now works in specific energies
and restores `state.e = floor * rho`.

### 2. the wall's admissibility bound read the energy wrongly

The bound compares the wall's power draw against the energy the near-wall
cell holds ABOVE the closure's floor. It computed that margin as:

```rust
let e_int = ctx.state.e[cell] - ke;    // J/m³ - J/kg
let margin_joules = (e_int - floor).max(0.0) * rho * volume;
```

Same mixing. Because `state.e` is about thirteen times the specific energy at
this density, the margin was inflated by roughly a factor of rho, and the
bound that was supposed to stop the wall stripping the cell barely engaged.
It is now:

```rust
let et = ctx.state.e[cell] / rho;
let e_int = et - ke;
```

### 3. the blowing-in film added kilograms to a density field

The film's injection accumulates a mass RATE into the state:

```rust
let dm = film_mass_rate[k] * dt;   // kg
state.rho[k] += dm;                // kg/m³ <- kg   WRONG
let de = film_energy_rate[k] * dt; // J
state.e[k] += de;                  // J/m³ <- J     WRONG
```

Both additions were missing the cell's volume, so a film injecting any
appreciable mass raised the density by hundreds of kilograms per cubic
metre in one step. That is precisely the observed signature: `rho_max`
exploding while every other field stayed physical. The species update had the
same problem in a third form, dividing by density where it needed the cell's
mass.

Now the amounts divide by the cell's volume, and the species enter as
`incoming species mass / cell mass`, which matches how the transport stage
treats the other species fluxes.

## Consequence for the wall bound

With the margins read correctly the bound finally does its job: a near-wall
cell that is approaching the floor has a shrinking margin, so the permitted
draw shrinks with it. A converged cell sitting well above the floor feels no
limit at all, which is what keeps the bound from distorting the answer.

The bound was also rewritten so that it distinguishes losses from gains. The
conduction term is signed along the face's axis while the radiant and film
terms are magnitudes of heat leaving the gas; turning all three into the same
outward sense lets the losing and gaining parts be separated. Heat INTO the
cell is never throttled, because giving a cell energy can never drive it
below the floor.

## Verification

- `cargo test -p nufor-core --lib`: 140 passed, 0 failed.
- `cargo clippy --workspace --all-targets -- -D warnings`: clean.
- The regression test `the_wall_bound_reads_the_cell_energy_as_specific_energy`
  was written to fail against the old code. It does: reverting the division
  reproduces the failure with "the withdrawal 275.16 exceeded the permitted
  quarter of the margin 120.96".

## Effect on earlier rungs

The guards were inert, so for runs that stayed healthy they changed nothing:
the comparison could not trip, so no cell was repaired and no cell was
damaged. Runs that DID cook now get a repair instead of silently continuing,
so their behaviour changes. Every rung whose gate depends on a cooled wall, a
strongly heated chamber or a combustion lane's residual history needs to be
re-run before its numbers are quoted again.

Cases affected by the change: the engine wall bracket's cooled legs, and any
case using `equations = "euler_axi"`, `combustion`, `cutcell` or the viscous
lane.

## A second audit: the species-enthalpy sign in the axisymmetric lane

The unit audit above fixed a dimensional error. A separate audit followed,
looking for the sign errors that hide behind a converging residual, and it
found one in the same family of code: the combustion lanes' coupling of
species diffusion to the energy balance.

### The convention

The conservation of total energy for a mixture with multispecies diffusion
carries the transported enthalpy on the diffusing mass:

    d(rho e)/dt + div(rho u e + p u) = -div(u.tau) - div(q) - div(sum_s h_s J_s)

Each combustion lane therefore stores an energy FLUX on the faces and applies
MINUS the divergence of what it stored. Writing `TE` for the stored flux, the
stored species term has to be `+ sum_s h_s J_s` for the update to come out
with the `-div(sum_s h_s J_s)` the law asks for. The Fourier term is stored
with its physical sign (`-k grad T`) for the same reason.

### What was found

Three lanes existed. Their face constructions turned out to be identical:
the same `js = -rho_f * d_f * dydx` species mass flux, the same
`qx = -k_f * dtdx`, the same `te = -(u.tau) + q +- sum(h_s J_s)` shape, and
the same `de = -(div TE)` update. Three of them stored `+ q_enth`:

    combustion_viscous.rs        + q_enth
    combustion_viscous_axi.rs    - q_enth     <- the odd one out
    cutcell_viscous.rs           + q_enth

The axisymmetric lane coupled the species diffusion to the energy with the
opposite sign from its two siblings, so wherever a composition gradient
existed the axial engine lane booked an energy change of the wrong
direction.

A fourth lane, `species_viscous.rs`, also stores `- q_enth`, but it applies
PLUS the divergence of its stored flux, so its arrangement is internally
consistent and it was left alone.

### How it was proved

The global sum of a divergence cannot reveal a sign error: the wrong sign
still sums to zero. A first attempt at a total-energy conservation test
passed against the defective code for exactly that reason, and a second
attempt that mixed in radial shear was measuring the axisymmetric
approximation rather than the sign. Both were discarded.

The decisive test compares the two lanes against each other. On a field that
varies only along the axis, every cylindrical term is identically zero: the
radial fluxes vanish, the radial velocity is zero, and the axial operator is
the planar operator. The two lanes therefore have to agree on the energy,
the momentum and the diffused species, cell for cell.

Against the defective code the energy changed by a factor of roughly five
between the lanes. With the sign corrected the two lanes agree to within
one part in ten billion across the whole row, and the momenta and the
species agree as well.

The test is `the_axisymmetric_lane_matches_the_planar_lane_on_an_axial_field`
in `crates/nufor-core/src/combustion_viscous_axi.rs`. It was fault-injected:
putting the old sign back makes it fail with a disagreement of 508 percent,
so it discriminates the defect rather than passing vacuously.

### Verification

- `cargo test -p nufor-core --lib`: 161 passed, 0 failed.
- `cargo clippy --workspace --all-targets -- -D warnings`: clean.
- Face-by-face comparison of the three combustion lanes shows identical
  expressions for the stresses, the Fourier flux and the species term.

### Effect on earlier rungs

The engine wall bracket runs on the axisymmetric lane, so every number
produced by the earlier bracket campaign is contaminated by this sign and has
to be recomputed. The bracket's gate is unaffected in its logic, but its
numbers are not quotable until the campaign is re-run.

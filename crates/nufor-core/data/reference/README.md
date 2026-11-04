# GRI-Mech 3.0 reference data

This directory holds the checked-in reference data behind the generated
GRI-Mech 3.0 tables in `crates/nufor-core/src/gri30_generated.rs` and the
parity fixture in `gri30_falloff_oracle.csv`.

## Provenance

- mechanism: CMU GRI-Mech 3.0
- source yaml: `/home/ebarleb/.local/lib/python3.12/site-packages/cantera/data/gri30.yaml`
- source sha256: `06650b1e0ee0012f6903d5328b1bb218cb6007d07f8ebe375d18f24811039345`
- cantera version: 3.2.0
- pyyaml version: 6.0.1
- species: 53
- reactions: 325, of which 29 are falloff reactions (26 troe, 3 lindemann)

The mechanism yaml is **not** vendored into this repository. Only derived
numeric tables and derived oracle values are checked in.

## Licensing

GRI-Mech 3.0 is distributed by MIT under the terms published by MIT
(`http://www.mit.edu/globalchemistry/`), which makes the mechanism freely
available subject to attribution. The license text itself is not
reproduced here. The checked-in files are derived numeric data (species
names and weights, reaction rate parameters, third-body efficiencies, and
oracle rate constants), not a copy of the licensed mechanism files.

## Unit conversions

The yaml is in cm-mol-s with activation energies in cal/mol. The tables
are in SI with m, kmol, s:

| quantity | yaml | table |
| --- | --- | --- |
| concentration | mol/cm3 | kmol/m3 |
| gas constant | 8314.46261815324 J/(kmol K) | 8314.46261815324 J/(kmol K) |
| activation energy | cal/mol | J/kmol, multiplied by 4184 |
| bimolecular rate constant | cm3/(mol s) | m3/(kmol s) |
| trimolecular rate constant | cm6/(mol2 s) | m6/(kmol2 s) |

Derivation of the pre-exponential factors, stated as the factor applied to
the yaml value:

- bimolecular, `cm3/(mol s) -> m3/(kmol s)`: one cubic centimetre is
  `1e-6 m3`, and one mole is `1e-3 kmol`, so the volume unit shrinks by
  `1e-6` while the amount unit grows by `1e-3`. Because the rate constant
  carries one inverse-volume factor per reacting molecule and the
  concentration moves from `mol/cm3` to `kmol/m3` by a factor of `1e3`,
  the combined factor is `1e-3`.
- trimolecular, `cm6/(mol2 s) -> m6/(kmol2 s)`: the six-length scale
  contributes `1e-12` and the squared amount contributes `1e-6`, so the
  combined factor is `1e-6`.
- unimolecular, `1/s -> 1/s`: unchanged.

These factors were verified numerically against Cantera's own converted
rate objects: for every reaction, the ratio of Cantera's converted
pre-exponential factor to the yaml value takes exactly one of `1e-3`,
`1e-6`, or `1`, selected by the reaction's molecularity.

Concentrations are quoted in kmol/m3, which is numerically equal to mol/L.
Note that `crates/nufor-core/src/combustion.rs` documents a different
historical convention for the Li 2004 mechanism, namely mol/m3 with
`E * 4184` in J/mol and a `1e-6` factor for bimolecular constants; the
GRI tables here therefore sit on the kmol basis and must not be mixed
with the Li tables without converting.

## Third-body convention

The bath concentration of a reaction is

```
[M] = sum over all species of efficiency_i * c_i
```

where `c_i` is the species molar concentration in kmol/m3, the
efficiency of a species named in the reaction's efficiency table is that
table's value, and the efficiency of a species the table omits is the
reaction's default efficiency.

There are two kinds of third body, and the default efficiency is what
distinguishes them:

- a generic `M` collider has a default efficiency of `1.0`, so every
  species contributes at least its own concentration, scaled by its
  table entry where the table names it.
- a named collider, spelled out in the reaction's equation, has a
  default efficiency of `0.0`, so only the named species contributes
  and every other species contributes nothing.

The named-collider reactions exist alongside the generic ones in this
mechanism and are different reactions with different rates. `H + O2 +
AR <=> HO2 + AR` uses argon alone, while `H + O2 + M <=> HO2 + M`
accepts every species with its own efficiency.

Two further rules apply:

- the bath multiplies both the forward and the reverse branch, so the
  same collider serves the reaction in either direction.
- the collider is not part of the stoichiometry. It cancels from both
  sides of the reaction, so the reactant and product lists exclude it
  even though the equation's text spells it out.

All four points were verified numerically against Cantera 3.2.0, whose
own `gas.third_body_concentrations` agrees to 2e-16 relative on the
generic reactions and reproduces the named-collider restriction
exactly.

## Troe formulation

For a troe falloff reaction, with `A`, `T3`, `T1`, `T2` taken from the
yaml as temperatures in kelvin:

```
Fc = (1 - A) exp(-T/T3) + A exp(-T/T1) + exp(-T2/T)
w0 = log10(max(Fc, 1e-300))
log10(pr) = log10(max(pr, 1e-300))
c = -0.4 - 0.67 w0
n = 0.75 - 1.27 w0
f1 = (log10(pr) + c) / (n - 0.14 (log10(pr) + c))
F  = 10^(w0 / (1 + f1^2))
```

The `exp(-T2/T)` term is additive and is present only when `T2` is
non-zero. The reduced pressure is dimensionless,

```
pr = [M] k0 / kinf
```

and the blended forward rate constant is

```
k = kinf (pr / (1 + pr)) F
```

A lindemann falloff reaction has `F = 1`. The Rust implementation lives
in `crates/nufor-core/src/troe.rs` and the mechanism wrapper in
`crates/nufor-core/src/gri30_mechanism.rs` drives every falloff reaction
through it.

### How the oracle removes the unit ambiguity

The fixture's broadening factor does not come from the physical reaction.
It comes from a deliberately constructed Cantera rate whose low- and
high-pressure branches are identical (`A = 1`, `b = 0`, `Ea = 0`). For
such a rate Cantera's reduced pressure is exactly the argument passed to
`falloff_function`, so the reference factor depends only on the Troe
formulation and not on any transcribed rate constant or on Cantera's
internal concentration scale. The fixture also carries the physical
reaction's reduced pressure, so the Rust side compares like with like.

The fixture additionally records Cantera's own forward rate constant and
the rate obtained by projecting Cantera's own `low_rate` and `high_rate`
through `k = kinf pr/(1+pr) F`; those two agree to 4.6e-15 relative and
serve as the secondary full-blend check.

## Fixture columns

`gri30_falloff_oracle.csv` holds one row per falloff reaction per
condition, covering 29 reactions x 6 temperatures x 3 pressures x 3
mixtures = 1566 rows plus a header:

| column | meaning |
| --- | --- |
| `reaction_index` | index into the mechanism's reaction table |
| `equation` | the reaction equation |
| `temperature_K` | temperature in kelvin |
| `pressure_Pa` | pressure in pascal |
| `mole_fraction_csv` | semicolon-separated `name:fraction` pairs |
| `efficiency_weighted_M_mol_per_L` | the bath from the convention above |
| `cantera_third_body_M_mol_per_L` | Cantera's own third-body concentration |
| `reduced_pressure_pr` | the dimensionless reduced pressure |
| `broadening_factor_cantera_clean` | Cantera's clean-oracle broadening factor |
| `k_fwd_cantera_SI` | Cantera's own forward rate constant |
| `k_fwd_projected_SI` | the projected blended forward rate constant |

Temperatures span 600 K to 3000 K, pressures span 0.1 atm to 10 atm, and
the mixtures include a hydrogen-steam blend (water at a non-default
collision efficiency of 6.0 in the affected reactions) and a dilute
radical blend whose argon is absent from most efficiency tables and so
contributes at the default efficiency of one.

## Regeneration

From the repository root:

```
env -u PYTHONPATH /usr/bin/python3.12 \
  crates/nufor-core/tools/generate_gri30_tables.py [yaml_path]
```

The script rewrites `crates/nufor-core/src/gri30_generated.rs` and
`crates/nufor-core/data/reference/gri30_falloff_oracle.csv`. It requires
Cantera 3.2.0 and PyYAML on the interpreter's path, and it prints the
sha256 of both files. Both outputs are deterministic: the same inputs
give byte-identical files, because the floats are emitted with Python's
shortest round-tripping representation and nothing carries a timestamp.

The thermodynamic coefficients come from a separate generator, which
rewrites `crates/nufor-core/src/gri30_thermo.rs`:

```
env -u PYTHONPATH /usr/bin/python3.12 tools/gen_gri30_thermo.py
```

The Cantera comparison fixtures are generated by
`tools/gen_gri30_netrates.py`, which writes the net-rate fixture, and by
`tools/gen_gri30_eqfixture.py`, which writes the equilibrium-constant
fixture.

## Scope

This directory holds reference data and provenance for the mechanism,
not a solver. The NASA-7 thermodynamic coefficients, the species names
and weights, the third-body efficiencies, the falloff parameters and the
equilibrium constants live in the generated Rust module and its
companion modules, not here.

The mechanism is selectable from a case file for a zero-dimensional
parcel. `cases/chemistry/gri30-stoich-ignition/case.toml` is marched by
`nufor ignite`, and the parcel's settled temperature and pressure match
Cantera's own constant-volume reactor to 0.0001 K and 0.01 Pa. The
comparison trace is generated by `tools/gen_ignition_reference.py` and
the figures by `cases/chemistry/ignition_figures.py`.

A spatial case cannot select the mechanism yet. The combustion closure
still assumes the nine-species hydrogen lane, so putting methane into a
flame on a grid is the outstanding work.

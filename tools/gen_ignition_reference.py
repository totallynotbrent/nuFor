#!/usr/bin/env python3
"""generates the cantera reference trace for an ignition case.

this is the oracle the ignition page quotes: cantera's own
constant-volume reactor, integrated with its own integrator, for the
same parcel the rust march burns. the reference is sampled densely
enough that interpolating it onto the rust march's sample grid costs
far less than the temperature difference being measured.

usage:

    env -u PYTHONPATH python3.12 gen_ignition_reference.py \
        <case.toml> <out.csv>

the case file names the mechanism yaml through the pinned source
tree, the parcel's initial temperature and pressure, the parcel's
composition (or stoichiometric fuel and air), and the horizon.
"""

import csv
import os
import sys

# cantera is installed for python3.12 on this machine.
sys.path.insert(0, "/home/ebarleb/.local/lib/python3.12/site-packages")

try:
    import tomllib
except ImportError:  # pragma: no cover - python below 3.11
    import tomli as tomllib  # type: ignore[no-redef]

import cantera as ct  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_MECH = os.path.join(
    REPO, "plans", "combustion_research", "sources", "gri30.yaml"
)

REPORTED = [
    "CH4",
    "O2",
    "H2O",
    "CO2",
    "CO",
    "H2",
    "O",
    "OH",
    "HO2",
    "H",
    "N2",
    "AR",
]


AIR_N2_OVER_O2 = 3.76
CH4_OXYGEN_PER_FUEL = 2.0


def stoich_air(gas, fuel):
    """the stoichiometric fuel and air mass fractions.

    one mole of the hydrocarbon fuel burns with two moles of oxygen,
    and dry air carries 3.76 moles of nitrogen per mole of oxygen.
    """
    mw = lambda name: gas.molecular_weights[gas.species_index(name)]
    masses = {
        fuel: mw(fuel),
        "O2": CH4_OXYGEN_PER_FUEL * mw("O2"),
        "N2": CH4_OXYGEN_PER_FUEL * AIR_N2_OVER_O2 * mw("N2"),
    }
    total = sum(masses.values())
    return {name: mass / total for name, mass in masses.items()}


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    case_path = argv[1]
    out_path = argv[2]

    with open(case_path, "rb") as handle:
        case = tomllib.load(handle)

    parcel = case["parcel"]
    mech = parcel.get("mechanism_yaml", DEFAULT_MECH)
    gas = ct.Solution(mech)

    if "composition" in parcel:
        composition = dict(parcel["composition"])
    else:
        composition = stoich_air(gas, parcel["fuel"])

    initial_t = float(parcel["initial_temperature_k"])
    initial_p = float(parcel["pressure_pa"])
    final_t = float(parcel["final_time_s"])

    gas.TPY = initial_t, initial_p, composition
    reactor = ct.IdealGasReactor(gas)
    network = ct.ReactorNet([reactor])

    # sample ten times finer than the march's grid so that the
    # interpolation onto that grid is negligible beside the measured
    # difference.
    stride = float(parcel["sample_interval_s"]) / 10.0

    rows = []
    previous = 0.0
    while network.time < final_t:
        network.step()
        time = network.time
        if time - previous >= stride or not rows:
            thermo = reactor.phase
            row = [time, thermo.T, thermo.P]
            for name in REPORTED:
                row.append(thermo.Y[thermo.species_index(name)])
            rows.append(row)
            previous = time

    # the file is written with unix line endings, so the csv module's
    # carriage-return default does not fight the repository's rule.
    with open(out_path, "w", newline="\n") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            ["time_s", "temperature_K", "pressure_Pa"]
            + [f"massfrac_{name}" for name in REPORTED]
        )
        writer.writerows(rows)

    thermo = reactor.phase
    print(f"wrote {out_path} ({len(rows)} rows)")
    print(f"final time {network.time:.6e} s")
    print(f"final temperature {thermo.T:.4f} K")
    print(f"final pressure {thermo.P:.4f} Pa")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

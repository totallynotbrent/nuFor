# the dissociating nitrogen cylinder gate: the far-field state, the shock
# standoff, the stagnation-line peak temperature and the atomic nitrogen
# read off the marched field, then compared against the reference values
# in reference.json with a percentage tolerance.
#
# the marched field itself is a run product and stays out of the
# repository. samples.py distils the slices this gate measures into
# samples.csv, which is committed, so the published numbers can be
# audited from a clean checkout.
#
# run from the case directory:
#   nufor run case.toml            # writes the vtk and the run log
#   python3 samples.py             # writes samples.csv
#   python3 gates.py               # writes gates.json
#
# every number comes from the marched field through the case's own mesh,
# body and initial condition. nothing is hardcoded and nothing is read
# from outside the repository.
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

CASE = os.path.join(HERE, "case.toml")
SAMPLES = os.path.join(HERE, "samples.csv")
REFERENCE = os.path.join(HERE, "reference.json")

# percentage bands each measured quantity must stay inside.
TOLERANCE = {
    "freestream_temperature_pct": 2.0,
    "freestream_pressure_pct": 2.0,
    "freestream_mach_pct": 2.0,
    "shock_standoff_pct": 3.0,
    "stagnation_peak_temperature_pct": 3.0,
    "wall_n_mole_fraction_pct": 5.0,
    "peak_n_mass_fraction_pct": 3.0,
}

# the closure's species order and their molar masses in kg per mol.
SPECIES = (("y_N2", 28.0134), ("y_O2", 31.9988), ("y_N", 15.9994),
           ("y_O", 15.9994), ("y_NO", 46.0055))


def read_case(path):
    """the case's own freestream, mesh and body."""
    text = open(path).read()
    vals = {}

    def grab(*keys):
        for key in keys:
            pat = r"(?m)^\s*" + key + r"\s*=\s*([-\d.eE+]+)"
            m = re.search(pat, text)
            if m:
                return float(m.group(1))
        return None

    vals["rho"] = grab("rho")
    vals["u"] = grab("u")
    vals["p"] = grab("p")
    vals["R"] = grab("gas_constant")
    vals["nx"] = grab("nx")
    vals["ny"] = grab("ny")
    vals["x0"] = grab("x0")
    vals["x1"] = grab("x1")
    vals["r"] = grab("r")
    vals["cx"] = grab("cx")
    if vals["cx"] is None:
        vals["cx"] = 0.0
    missing = [k for k in vals if vals[k] is None]
    if missing:
        raise ValueError(f"{path}: missing {missing}")
    return vals


def read_samples(path):
    """the distiled field slices: the freestream slab, the stagnation
    line and the atomic nitrogen map of the whole field."""
    sections = {"inlet": [], "stagnation": [], "field_n": []}
    current = None
    header = None
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if not row:
                continue
            if row[0].startswith("#"):
                if "section" in row[0] and len(row) > 1:
                    current = row[1]
                continue
            if row[0] == "i" and current in ("inlet", "stagnation"):
                header = row
                continue
            if row[0] == "j" and current == "field_n":
                header = row[1:]
                continue
            if current is None or header is None:
                continue
            if current == "field_n":
                sections[current].append(
                    [float(v) for v in row[1].split(",")] if len(row) == 2 else
                    [float(v) for v in row[1:]]
                )
            else:
                sections[current].append(dict(zip(header, (float(v) for v in row))))
    for name, rows in sections.items():
        if not rows:
            raise ValueError(f"{path}: the {name} section is empty")
    return sections


def median(values):
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def mole_fraction(cell):
    """the mole fraction of atomic nitrogen in the mixture of a cell."""
    inv = 0.0
    for name, mass in SPECIES:
        y = cell.get(name)
        if y is None or y <= 0.0:
            continue
        inv += y / mass
    if inv <= 0.0:
        return 0.0
    return (cell["y_N"] / 15.9994) / inv


def measure(sections, case):
    """the published quantities, read off the distiled field."""
    inlet = sections["inlet"]
    stagnation = sections["stagnation"]

    def column(rows, name):
        return [row[name] for row in rows]

    freestream_temperature = median(column(inlet, "temperature"))
    freestream_pressure = median(column(inlet, "pressure"))
    freestream_mach = median(column(inlet, "mach"))

    pitch = (case["x1"] - case["x0"]) / case["nx"]

    def centre(i):
        return case["x0"] + (i + 0.5) * pitch

    # the stagnation line is the centreline row, and the body occupies
    # the chord around the sphere centre on it. the flow comes from the
    # west, so the leading edge of that chord is the nose.
    nose_x = case["cx"] - case["r"]
    nose_index = None
    for row in stagnation:
        if centre(int(row["i"])) >= nose_x:
            nose_index = int(row["i"])
            break
    if nose_index is None:
        raise ValueError("no nose cell on the stagnation line")

    # the shock: walking back from the nose, the pressured pocket runs
    # to the shock, and the first cell back at freestream level marks it.
    freestream_level = 1.05 * freestream_pressure
    shock_x = None
    for row in reversed(stagnation):
        i = int(row["i"])
        if i >= nose_index:
            continue
        if row["pressure"] < freestream_level:
            shock_x = centre(i + 1)
            break
    if shock_x is None:
        raise ValueError("no shock found on the stagnation line")
    shock_standoff = nose_x - shock_x

    # the peak temperature of the shock layer on the stagnation line.
    stagnation_peak_temperature = max(
        row["temperature"] for row in stagnation
        if abs(centre(int(row["i"])) - case["cx"]) >= case["r"]
    )

    # the last fluid cell ahead of the nose carries the near-wall
    # composition.
    wall_row = None
    for row in stagnation:
        if int(row["i"]) == nose_index - 1:
            wall_row = row
            break
    if wall_row is None:
        raise ValueError("no wall cell on the stagnation line")
    wall_n_mole_fraction = mole_fraction(wall_row)

    peak_n_mass_fraction = max(max(row) for row in sections["field_n"])

    return {
        "freestream_temperature_k": freestream_temperature,
        "freestream_pressure_pa": freestream_pressure,
        "freestream_mach": freestream_mach,
        "shock_standoff_m": shock_standoff,
        "stagnation_peak_temperature_k": stagnation_peak_temperature,
        "wall_n_mole_fraction": wall_n_mole_fraction,
        "peak_n_mass_fraction": peak_n_mass_fraction,
    }


def reference_key(key):
    return {
        "freestream_temperature_pct": "freestream_temperature_k",
        "freestream_pressure_pct": "freestream_pressure_pa",
        "freestream_mach_pct": "freestream_mach",
        "shock_standoff_pct": "shock_standoff_m",
        "stagnation_peak_temperature_pct": "stagnation_peak_temperature_k",
        "wall_n_mole_fraction_pct": "wall_n_mole_fraction",
        "peak_n_mass_fraction_pct": "peak_n_mass_fraction",
    }[key]


def main():
    if not os.path.exists(SAMPLES):
        print(f"MISS: no {os.path.basename(SAMPLES)}, run samples.py after the case")
        return 1
    case = read_case(CASE)
    sections = read_samples(SAMPLES)
    got = measure(sections, case)
    entry = {"tag": "park5", **got}
    failures = []

    # the far field must be the state the case declares. this is an
    # absolute invariant, not a drift check against the last run.
    declared_temperature = case["p"] / (case["rho"] * case["R"])
    if abs(got["freestream_temperature_k"] - declared_temperature) / declared_temperature > 0.02:
        failures.append(
            f"far field temperature {got['freestream_temperature_k']:.1f} K against the "
            f"declared {declared_temperature:.1f} K"
        )

    if os.path.exists(REFERENCE):
        reference = json.load(open(REFERENCE))[0]
        for key, tol in TOLERANCE.items():
            ref_key = reference_key(key)
            if ref_key not in reference:
                continue
            base = reference[ref_key]
            drift = abs(got[ref_key] - base) / abs(base) * 100.0 if base else 0.0
            entry[f"{key}_drift_pct"] = drift
            if drift > tol:
                failures.append(
                    f"{ref_key} {got[ref_key]:.6g} vs reference {base:.6g} "
                    f"({drift:.2f}% > {tol}%)"
                )

    with open(os.path.join(HERE, "gates.json"), "w") as fh:
        json.dump([entry], fh, indent=1)
        fh.write("\n")

    print(
        f"park5: far field {got['freestream_temperature_k']:.1f} K  "
        f"{got['freestream_pressure_pa']:.3f} Pa  Mach {got['freestream_mach']:.2f}  "
        f"standoff {got['shock_standoff_m']:.3f} m  peak T "
        f"{got['stagnation_peak_temperature_k']:.0f} K  wall X_N "
        f"{got['wall_n_mole_fraction']:.3f}  peak y_N {got['peak_n_mass_fraction']:.3f}"
    )
    if failures:
        print("")
        for f in failures:
            print(f"FAIL {f}")
        return 1
    print("")
    print("park5 gates: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""the chamber-to-throat gate.

the settled state of the chamber march ships with the case as a small
artifact, so this gate checks three things:

- the contracts hold: the throat chokes, the sonic pressure and
  temperature sit inside their tolerance around CEA, the mass flow is
  conserved, the chamber runs subsonic and the exit runs supersonic.
- the gas agrees with cea: the marched gas's own isentrope and cea's
  equilibrium throat have to sit close together, or the comparison is
  not about the gas.
- the march agrees with its gas: the choked pressure has to land near
  the isentrope the crate's heats give, within the band the marching
  error is allowed.

the page is read and compared with the recorded state, so a number
cannot drift away from the march that produced it without the gate
noticing.

the walked field itself stays out of the repository, so the gate also
checks the recorded state against the field whenever the field is on
disk (after a local march), and says so when it cannot.

usage:
    python3 gates.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

from chamber_fields import (  # noqa: E402
    read_case,
    read_vtk,
    sonic_point,
)

SETTLED = os.path.join(HERE, "settled")
RECORD = os.path.join(SETTLED, "settled.json")
RESULT = os.path.join(SETTLED, "thermo_result.csv")
CURVE = os.path.join(SETTLED, "isentrope.csv")
FIELD = os.path.join(SETTLED, "vtk-final.vtk")
PAGE = os.path.join(ROOT, "docs", "numerics", "2d", "chamber-throat.md")


def recorded():
    """the settled state, as the march recorded it."""
    with open(RECORD) as handle:
        state = json.load(handle)
    values = {
        "pressure_mpa": state["sonic"]["p_pa_mean"] / 1e6,
        "pressure_band_mpa": state["sonic"]["p_pa_sigma"] / 1e6,
        "temperature_k": state["sonic"]["T_k_mean"],
        "temperature_band_k": state["sonic"]["T_k_sigma"],
        "mass_flow_kg_s": state["sonic"]["m_dot_kg_s_mean"],
        "mass_flow_band_kg_s": state["sonic"]["m_dot_kg_s_sigma"],
        "chamber_pressure_mpa": state["chamber"]["p_pa_mean"] / 1e6,
        "chamber_pressure_band_mpa": state["chamber"]["p_pa_sigma"] / 1e6,
        "chamber_temperature_k": state["chamber"]["T_k_mean"],
        "chamber_temperature_band_k": state["chamber"]["T_k_sigma"],
    }
    profile = state["profile"]
    values["exit_mach"] = profile[-1]["flux_mach"]
    values["chamber_mach"] = profile[len(profile) // 8]["mixed_mach"]
    return values, state["cea"]


def pairs(path):
    """a two-column csv of key,value pairs."""
    out = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            key, _, value = line.partition(",")
            out[key.strip()] = value.strip()
    return out


def curve_rows():
    """the isentrope curve the rust audit wrote."""
    with open(CURVE) as handle:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(row for row in handle if not row.startswith("#"))
        ]


def contracts_from(values, cea, result):
    """the contract margins, as (name, margin, bound, side) tuples."""
    isentrope = float(result["isentrope_sonic_p_pa"])
    gas_vs_cea = abs(100.0 * (cea["throat_p_pa"] / isentrope - 1.0))
    march_vs_gas = abs(100.0 * (values["pressure_mpa"] * 1e6 / isentrope - 1.0))
    return (
        ("the sonic temperature", abs(100.0 * (values["temperature_k"] / cea["throat_T_k"] - 1.0)), 5.0, "at most"),
        ("the sonic pressure", abs(100.0 * (values["pressure_mpa"] * 1e6 / cea["throat_p_pa"] - 1.0)), 12.0, "at most"),
        ("the mass flow is conserved", abs(100.0 * (values["mass_flow_kg_s"] - 3.80) / 3.80), 8.0, "at most"),
        ("the chamber runs subsonic", values["chamber_mach"], 0.1, "below"),
        ("the exit runs supersonic", values["exit_mach"], 1.0, "above"),
        ("the gas agrees with cea", gas_vs_cea, 1.0, "at most"),
        ("the march agrees with its gas", march_vs_gas, 10.0, "at most"),
    )


# the page's measured table, label to (value, band, tolerance)
PUBLISHED = {
    "sonic pressure": ("pressure_mpa", "pressure_band_mpa", 0.002),
    "sonic temperature": ("temperature_k", "temperature_band_k", 0.5),
    "throat mass flow": ("mass_flow_kg_s", "mass_flow_band_kg_s", 0.002),
    "chamber pressure": ("chamber_pressure_mpa", "chamber_pressure_band_mpa", 0.002),
    "chamber temperature": ("chamber_temperature_k", "chamber_temperature_band_k", 0.5),
    "exit flux mach": ("exit_mach", None, 0.002),
}

# the page's split table, label to (value, tolerance)
PUBLISHED_SPLIT = {
    "the marched gas's own isentrope": ("isentrope_sonic_p_pa", 0.002),
    "cea's equilibrium throat": ("cea_throat_p_pa", 0.002),
    "the march": ("measured_sonic_p_pa", 0.002),
}


def published_rows():
    """the page's tables, as {label: [cells]}."""
    text = open(PAGE).read()
    rows = {}
    for line in text.splitlines():
        if not line.startswith("| ") or line.startswith("|--"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 2:
            continue
        label = cells[0].lower()
        rows.setdefault(label, []).append(cells)
    return rows


def first_number(cell):
    try:
        return float(cell.split()[0])
    except (ValueError, IndexError):
        return None


def main():
    failures = []
    values, cea = recorded()
    result = pairs(RESULT)

    # the march must have found a sonic point at all.
    history = list(csv.DictReader(open(os.path.join(SETTLED, "history.csv"))))
    pressures = [float(row["sonic_p_pa"]) for row in history]
    if not pressures:
        print("FAIL the record holds no sonic point")
        failures.append("the sonic point")
    elif pressures[0] == pressures[-1]:
        print("FAIL the sonic point never moved, the march did not settle")
        failures.append("the settling")

    # the isentrope curve has to reach the sonic point and lose pressure
    # on the way there.
    curve = curve_rows()
    if not curve:
        print("FAIL the isentrope curve is missing")
        failures.append("the isentrope curve")
    else:
        loss = curve[-1]["deviation_percent"]
        if not -12.0 <= loss <= 0.0:
            print(f"FAIL the section's loss sits at {loss:.2f} percent")
            failures.append("the section's loss")
        else:
            print(f"PASS the section's loss sits at {loss:.2f} percent")

    if os.path.exists(FIELD):
        case = read_case(os.path.join(HERE, "case.toml"))
        fields, nx, ny, origin, spacing = read_vtk(FIELD)
        sonic = sonic_point(case, fields, nx, ny, origin, spacing)
        if sonic is None:
            print("FAIL the field on disk carries no sonic point")
            failures.append("the field's sonic point")
        elif abs(sonic["p"] - values["pressure_mpa"] * 1e6) > 0.05e6:
            print(
                f"FAIL the field's sonic pressure {sonic['p']/1e6:.4f} MPa disagrees "
                f"with the record {values['pressure_mpa']:.4f} MPa"
            )
            failures.append("the recorded state against the field")
        else:
            print("PASS the recorded state agrees with the field")
    else:
        print("SKIP the field is not on disk, so the recorded state stands unchecked against it")

    for name, margin, bound, side in contracts_from(values, cea, result):
        if side == "at most":
            ok = margin <= bound
        elif side == "below":
            ok = margin < bound
        else:
            ok = margin > bound
        print(("PASS " if ok else "FAIL ") + f"{name}: {margin:.4f} against {bound}")
        if not ok:
            failures.append(name)

    rows = published_rows()
    missing = sorted(set(PUBLISHED) - set(rows))
    if missing:
        failures.append("the page's measured table is missing rows: " + ", ".join(missing))
        print("FAIL the page's measured table is missing rows: " + ", ".join(missing))
    for label, (key, band_key, tol) in PUBLISHED.items():
        if label not in rows:
            continue
        cells = rows[label][0]
        got = first_number(cells[1])
        want = values[key]
        ok = got is not None and abs(got - want) <= tol
        print(("PASS " if ok else "FAIL ") + f"the page's {label}: {got} against the recorded {want:.4f}")
        if not ok:
            failures.append(f"the page's {label}")
        if band_key and len(cells) > 2:
            got_band = first_number(cells[2])
            want_band = values[band_key]
            ok = got_band is not None and abs(got_band - want_band) <= tol
            print(("PASS " if ok else "FAIL ")
                  + f"the page's {label} band: {got_band} against the recorded {want_band:.4f}")
            if not ok:
                failures.append(f"the page's {label} band")

    split_missing = sorted(set(PUBLISHED_SPLIT) - set(rows))
    if split_missing:
        failures.append("the page's split table is missing rows: " + ", ".join(split_missing))
        print("FAIL the page's split table is missing rows: " + ", ".join(split_missing))
    split_values = {
        "isentrope_sonic_p_pa": float(result["isentrope_sonic_p_pa"]) / 1e6,
        "cea_throat_p_pa": float(cea["throat_p_pa"]) / 1e6,
        "measured_sonic_p_pa": values["pressure_mpa"],
    }
    for label, (key, tol) in PUBLISHED_SPLIT.items():
        if label not in rows:
            continue
        got = first_number(rows[label][0][1])
        want = split_values[key]
        ok = got is not None and abs(got - want) <= tol
        print(("PASS " if ok else "FAIL ") + f"the page's {label}: {got} against the recorded {want:.4f}")
        if not ok:
            failures.append(f"the page's {label}")

    if failures:
        print("FAILED: " + "; ".join(failures))
        return 1
    print("passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

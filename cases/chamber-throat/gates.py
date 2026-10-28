"""the chamber-to-throat gate.

the settled state of the chamber march ships with the case as a small
artifact, so this gate checks two things:

- the contracts hold: the throat chokes, the throat's pressure and
  temperature sit inside their tolerance around CEA, the mass flow is
  conserved, the chamber runs subsonic and the exit runs supersonic.
- the page publishes those numbers. the page's tables are read and
  compared with the recorded state, so a number cannot drift away from
  the march that produced it without the gate noticing.

the walked field itself stays out of the repository, so the gate also
checks the recorded state against the field whenever the field is on
disk (after a local march), and says so when it cannot.

usage:
    python3 gates.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)

from chamber_fields import (  # noqa: E402
    read_case,
    read_vtk,
    columns,
    throat_column,
)

SETTLED = os.path.join(HERE, "settled")
RECORD = os.path.join(SETTLED, "settled.json")
FIELD = os.path.join(SETTLED, "vtk-final.vtk")
PAGE = os.path.join(ROOT, "docs", "numerics", "2d", "chamber-throat.md")


def recorded():
    """the settled state, as the march recorded it."""
    with open(RECORD) as handle:
        state = json.load(handle)
    return {
        "mach": state["throat_mach"],
        "pressure_mpa": state["throat_p_pa"] / 1e6,
        "temperature_k": state["throat_T_k"],
        "mass_flow_kg_s": state["throat_mdot_kg_s"],
        "chamber_pressure_mpa": state["chamber_p_pa"] / 1e6,
        "chamber_temperature_k": state["chamber_T_k"],
        "exit_mach": state["exit_mach"],
        "chamber_mach": state["chamber_mach"],
        "throat_mdot": state["throat_mdot_kg_s"],
        "chamber_mdot": state["chamber_mdot_kg_s"],
    }, state["cea"]


def from_field(state):
    """the exit state, recomputed from the field when the field is on disk."""
    case = read_case(os.path.join(HERE, "case.toml"))
    fields, nx, ny, origin, spacing = read_vtk(FIELD)
    cols = columns(case, fields, nx, ny, origin, spacing)
    _, _ = throat_column(case, cols)
    exit_col = sorted(cols)[-1]
    values = dict(state)
    values["exit_mach"] = cols[exit_col]["mach"]
    return values


def contracts_from(values, cea):
    """the contract margins, as (name, margin, bound, above) tuples."""
    return (
        ("the throat chokes", abs(values["mach"] - 1.0), 0.08, "at most"),
        ("the throat's temperature", abs(100.0 * (values["temperature_k"] / cea["throat_T_k"] - 1.0)), 5.0, "at most"),
        ("the throat's pressure", abs(100.0 * (values["pressure_mpa"] * 1e6 / cea["throat_p_pa"] - 1.0)), 12.0, "at most"),
        ("the mass flow is conserved", abs(100.0 * (values["chamber_mdot"] - values["throat_mdot"]) / values["throat_mdot"]), 8.0, "at most"),
        ("the chamber runs subsonic", values["chamber_mach"], 0.1, "below"),
        ("the exit runs supersonic", values["exit_mach"], 1.0, "above"),
    )


PUBLISHED = {
    "mach": ("mach", 0.002),
    "pressure": ("pressure_mpa", 0.002),
    "temperature": ("temperature_k", 0.5),
    "mass flow": ("mass_flow_kg_s", 0.002),
    "chamber pressure": ("chamber_pressure_mpa", 0.002),
    "chamber temperature": ("chamber_temperature_k", 0.5),
    "exit mach": ("exit_mach", 0.002),
}


def published_rows():
    """the page's measured table, as {label: value}."""
    text = open(PAGE).read()
    rows = {}
    for line in text.splitlines():
        if not line.startswith("| ") or line.startswith("|--"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 2:
            continue
        label = cells[0].lower()
        if label in PUBLISHED:
            try:
                rows[label] = float(cells[1].split()[0])
            except ValueError:
                pass
    return rows


def main():
    failures = []
    values, cea = recorded()

    if os.path.exists(FIELD):
        field_values = from_field(values)
        drifted = [
            key for key, tol in (("mach", 0.002), ("exit_mach", 0.002))
            if abs(field_values[key] - values[key]) > tol
        ]
        if drifted:
            print("FAIL the recorded state disagrees with the field on: " + ", ".join(drifted))
            failures.append("the recorded state against the field")
        else:
            print("PASS the recorded state agrees with the field")
    else:
        print("SKIP the field is not on disk, so the recorded state stands unchecked against it")

    for name, margin, bound, side in contracts_from(values, cea):
        if side == "at most":
            ok = margin <= bound
        elif side == "below":
            ok = margin < bound
        else:
            ok = margin > bound
        print(("PASS " if ok else "FAIL ") + f"{name}: {margin:.4f} against {bound}")
        if not ok:
            failures.append(name)

    published = published_rows()
    missing = sorted(set(PUBLISHED) - set(published))
    if missing:
        failures.append("the page's measured table is missing rows: " + ", ".join(missing))
        print("FAIL the page's measured table is missing rows: " + ", ".join(missing))
    for label, (key, tol) in PUBLISHED.items():
        if label not in published:
            continue
        want = values[key]
        got = published[label]
        ok = abs(got - want) <= tol
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

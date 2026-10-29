"""compares the chamber march across grid resolutions.

each record directory is one distilled level: the distiller's history,
profile and thermo inputs, plus the rust audit's thermo result. the audit
runs once per record with NUFOR_CHAMBER_RECORD pointing at it, so every
level carries the same isentrope walk.

the comparison answers one question: does the marcher's own deficit
shrink as the grid tightens, and does the gas-versus-cea term stay small
so the trend is about the meshing and not about the reference.

usage:
    python3 resolution.py <record dir> [<record dir> ...]
    python3 resolution.py --write <summary csv> <record dir> [...]
    python3 resolution.py --summary <summary csv>
"""
import csv
import json
import os
import sys

# the largest deficit the marcher is allowed to owe its own gas
DEFICIT_BOUND_PERCENT = 10.0
# the gas against cea has to stay inside this on every level
GAS_VS_CEA_BOUND_PERCENT = 1.0

SUMMARY_COLUMNS = (
    "grid",
    "nx",
    "ny",
    "sonic_p_mpa",
    "isentrope_p_mpa",
    "deficit_percent",
    "gas_vs_cea_percent",
    "sonic_x_mm",
    "mass_flow_kg_s",
)


def read_pairs(path):
    """reads a two-column csv of key,value pairs."""
    out = {}
    with open(path, newline="") as handle:
        for row in csv.reader(handle):
            if len(row) < 2:
                continue
            key = row[0].strip()
            if not key or key.startswith("#"):
                continue
            out[key] = row[1].strip()
    return out


def read_level(record_dir):
    """one resolution level's numbers, read from the record itself."""
    record = json.load(open(os.path.join(record_dir, "settled.json")))
    result = read_pairs(os.path.join(record_dir, "thermo_result.csv"))
    inputs = read_pairs(os.path.join(record_dir, "thermo_inputs.csv"))

    deficit = float(result["measured_minus_isentrope_percent"])
    isentrope_p = float(result["isentrope_sonic_p_pa"])
    cea_p = float(inputs["cea_throat_p_pa"])
    measured_p = record["sonic"]["p_pa_mean"]
    gas_gap = 100.0 * (cea_p / isentrope_p - 1.0)

    nx = record["grid"][0]
    ny = record["grid"][1]
    return {
        "grid": f"{nx}x{ny}",
        "nx": nx,
        "ny": ny,
        "cells": nx * ny,
        "sonic_p_mpa": measured_p / 1.0e6,
        "isentrope_p_mpa": isentrope_p / 1.0e6,
        "deficit_percent": deficit,
        "gas_vs_cea_percent": gas_gap,
        "sonic_x_mm": record["sonic"]["x_m"] * 1000.0,
        "mass_flow_kg_s": record["sonic"]["m_dot_kg_s_mean"],
    }


def read_summary(path):
    """the committed summary: one row per grid level."""
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            grid = row["grid"].strip()
            nx, _, ny = grid.partition("x")
            rows.append({
                "grid": grid,
                "nx": int(nx),
                "ny": int(ny),
                "cells": int(nx) * int(ny),
                "sonic_p_mpa": float(row["sonic_p_mpa"]),
                "isentrope_p_mpa": float(row["isentrope_p_mpa"]),
                "deficit_percent": float(row["deficit_percent"]),
                "gas_vs_cea_percent": float(row["gas_vs_cea_percent"]),
                "sonic_x_mm": float(row["sonic_x_mm"]),
                "mass_flow_kg_s": float(row["mass_flow_kg_s"]),
            })
    return rows


def write_summary(path, levels):
    """writes the levels as the summary csv."""
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(SUMMARY_COLUMNS)
        for level in sorted(levels, key=lambda item: item["cells"]):
            writer.writerow([
                level["grid"],
                level["nx"],
                level["ny"],
                f"{level['sonic_p_mpa']:.6f}",
                f"{level['isentrope_p_mpa']:.6f}",
                f"{level['deficit_percent']:.6f}",
                f"{level['gas_vs_cea_percent']:.6f}",
                f"{level['sonic_x_mm']:.4f}",
                f"{level['mass_flow_kg_s']:.6f}",
            ])


def check_levels(levels):
    """the resolution contracts, as a list of failure strings."""
    failures = []
    previous = None
    for level in sorted(levels, key=lambda item: item["cells"]):
        if previous is not None:
            change = abs(level["deficit_percent"]) - abs(previous["deficit_percent"])
            if change > 0.0:
                failures.append(
                    f"{level['grid']} does not beat {previous['grid']} on the marcher's deficit"
                )
        if abs(level["deficit_percent"]) > DEFICIT_BOUND_PERCENT:
            failures.append(
                f"{level['grid']} leaves the marcher owing {level['deficit_percent']:.2f} percent"
            )
        if abs(level["gas_vs_cea_percent"]) > GAS_VS_CEA_BOUND_PERCENT:
            failures.append(
                f"{level['grid']} disagrees with cea by {level['gas_vs_cea_percent']:.2f} percent"
            )
        previous = level

    ordered = sorted(levels, key=lambda item: item["cells"])
    if len(ordered) >= 2:
        first = ordered[0]["deficit_percent"]
        last = ordered[-1]["deficit_percent"]
        if abs(last) >= abs(first):
            failures.append("the marcher's own deficit did not shrink with the grid")
    return failures


def report(levels):
    """prints the table."""
    print(
        f"{'grid':>12} {'sonic MPa':>10} {'isentrope':>10} {'deficit':>9} "
        f"{'gas vs cea':>11} {'delta':>8}"
    )
    previous = None
    for level in sorted(levels, key=lambda item: item["cells"]):
        delta = ""
        if previous is not None:
            change = abs(level["deficit_percent"]) - abs(previous["deficit_percent"])
            delta = f"{change:+.2f} pp"
        print(
            f"{level['grid']:>12} {level['sonic_p_mpa']:>10.4f} "
            f"{level['isentrope_p_mpa']:>10.4f} {level['deficit_percent']:>+8.2f}% "
            f"{level['gas_vs_cea_percent']:>+10.2f}% {delta:>8}"
        )
        previous = level
    ordered = sorted(levels, key=lambda item: item["cells"])
    if len(ordered) >= 2:
        print(
            f"the marcher's own deficit {ordered[0]['deficit_percent']:+.2f} percent -> "
            f"{ordered[-1]['deficit_percent']:+.2f} percent"
        )


def main():
    argv = sys.argv[1:]
    if not argv:
        raise SystemExit("give at least one distilled record directory")

    if argv[0] == "--summary":
        if len(argv) != 2:
            raise SystemExit("--summary takes one summary csv")
        levels = read_summary(argv[1])
    else:
        write = None
        if argv[0] == "--write":
            if len(argv) < 3:
                raise SystemExit("--write takes a summary csv and at least one record directory")
            write = argv[1]
            argv = argv[2:]
        levels = []
        for record_dir in argv:
            try:
                levels.append(read_level(record_dir))
            except FileNotFoundError as err:
                raise SystemExit(
                    f"{record_dir} is missing {err.filename}: run the distiller and the "
                    "rust audit for that record first"
                )
        if write:
            write_summary(write, levels)

    report(levels)
    failures = check_levels(levels)
    for failure in failures:
        print(f"FAIL {failure}")
    if failures:
        raise SystemExit(1)
    print("passed")


if __name__ == "__main__":
    main()

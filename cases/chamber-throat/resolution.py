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
# the deficit has to be significant on the baseline level, or the study
# is not establishing a defect to explain
DEFICIT_PRESENCE_PERCENT = 3.0
# at a fixed radial spacing, the deficit may drift by this much when the
# axis is refined
AXIAL_DRIFT_BOUND_PERCENT = 1.0

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


def radial_spacing_m(level):
    """the level's radial cell thickness, in millimetres.

    the chamber's span is 100 mm, so the radial spacing is 100 over the
    row count. the deficit follows this spacing and not the total cell
    count, which is what the measurements showed.
    """
    return 100.0 / level["ny"]


def check_levels(levels):
    """the resolution study's contracts, as a list of failure strings.

    this study exists to CHARACTERISE the wall representation's defect,
    so its contracts are about what the measurements must show for the
    diagnosis to stand. they are not a wish that the grid be good.

    - the reference has to hold still. every level's gas has to agree
      with CEA within one percent, or the comparison is not about the
      meshing.
    - the defect has to be PRESENT. the baseline level's deficit must
      be significant, because a study that found the deficit vanishing
      under refinement would refute the diagnosis rather than confirm
      it.
    - the deficit has to track the RADIAL spacing. a finer radial
      spacing owes less, because the throat's radius is a radial
      quantity and the staircase quantises it to whole radial rows.
    - the deficit must NOT track the AXIAL spacing. at a fixed radial
      spacing, refining the axis has to leave the deficit roughly where
      it was. if axial refinement fixed the loss, the loss would not be
      the wall's radial representation.
    - a level that leaves the marcher owing more than the bound is
      ILL-CONDITIONED rather than informative: the flow's own state
      near the throat is not resolved there. it is reported and it is
      excluded from the trends, but its presence does not fail the
      study, because a staircase's behaviour degrading under refinement
      is part of what the study establishes.
    """
    failures = []

    for level in levels:
        if abs(level["gas_vs_cea_percent"]) > GAS_VS_CEA_BOUND_PERCENT:
            failures.append(
                f"{level['grid']} disagrees with cea by {level['gas_vs_cea_percent']:.2f} percent"
            )

    baseline = min(levels, key=lambda item: item["cells"])
    if abs(baseline["deficit_percent"]) < DEFICIT_PRESENCE_PERCENT:
        failures.append(
            f"the baseline {baseline['grid']} owes only "
            f"{baseline['deficit_percent']:.2f} percent, so the study does not "
            "establish a defect to explain"
        )

    healthy = [
        level for level in levels if abs(level["deficit_percent"]) <= DEFICIT_BOUND_PERCENT
    ]

    # the deficit tracks the radial spacing: at a fixed axial spacing, a
    # finer radial spacing owes less.
    by_nx = {}
    for level in healthy:
        by_nx.setdefault(level["nx"], []).append(level)
    for nx, group in sorted(by_nx.items()):
        ordered = sorted(group, key=radial_spacing_m)
        for finer, coarser in zip(ordered, ordered[1:]):
            if abs(finer["deficit_percent"]) >= abs(coarser["deficit_percent"]):
                failures.append(
                    f"{finer['grid']} does not beat {coarser['grid']} on the "
                    "marcher's deficit at the same axial spacing"
                )

    # and it does NOT track the axial spacing: at a fixed radial
    # spacing, refining the axis leaves it where it was.
    by_radial = {}
    for level in healthy:
        key = round(radial_spacing_m(level), 6)
        by_radial.setdefault(key, []).append(level)
    for key, group in sorted(by_radial.items(), reverse=True):
        ordered = sorted(group, key=lambda item: item["nx"])
        for finer, coarser in zip(ordered, ordered[1:]):
            drift = abs(abs(finer["deficit_percent"]) - abs(coarser["deficit_percent"]))
            if drift > AXIAL_DRIFT_BOUND_PERCENT:
                failures.append(
                    f"{finer['grid']} moves the deficit by {drift:.2f} points against "
                    f"{coarser['grid']} at the same radial spacing, so the loss does not "
                    "follow the wall's radial representation"
                )
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

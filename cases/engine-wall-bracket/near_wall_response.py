#!/usr/bin/env python3
"""the near-wall gas temperature: does the wall reach the gas it touches?

WHY THIS EXISTS

the wall's heat flux has to reach the gas on the other side of it. the
visible symptoms of a wall that does not reach its gas are a believable
reported duty and a gas that does not care, so the interesting
measurement is the gas's response, not the wall's claim.

this module measures the temperature in the cells ADJACENT TO THE WALL
and compares the wall's legs against the adiabatic leg. the near-wall
gas is where the wall's influence must show up first and strongest.

WHAT IT CAUGHT

the viscous axisymmetric lane's radial divergence read the first column's
radial faces for every column. the wall's heat flux was therefore spent
on the gas above the inlet, and the gas under the wall saw nothing. the
reported duty stayed believable while the gas did not respond to it.

WHY THERE IS NO FLOW PERIOD HERE

an earlier version of this module averaged the near-wall gas over a
presumed flow period. measurement refuted the premise: the saved fields
do not repeat, so an average over a guessed period carries the wandering
phase along with it and washes out the very difference it was meant to
expose. the comparison is now made the same way the gate makes it: the
near-wall temperature of each leg is differenced against the adiabatic
leg at the SAME saved field, over the window the gate's settlement
criterion chose. the pairing cancels the chamber's own motion, which
is common to every leg and owned by none.

USAGE
    python3 near_wall_response.py <run_root>
"""

from __future__ import annotations

import glob
import importlib.util
import json
import os
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = HERE.parent.parent / "cases"


def load_field_reader():
    """the chamber-throat field reader, which reads the marched vtks."""
    candidate = os.path.join(str(CASES), "chamber-throat")
    sys.path.insert(0, candidate)
    import chamber_fields

    return chamber_fields


def near_wall_cells(mask, nx: int, ny: int) -> list[int]:
    """the gas cells that touch a solid cell.

    a gas cell qualifies when it shares a face with the solid: directly
    above, below or beside it. those are the cells the wall's flux acts
    on, and they are where the wall's effect must appear.
    """
    near = []
    for j in range(ny):
        for i in range(nx):
            k = j * nx + i
            if mask[k]:
                continue
            neighbours = ((j + 1, i), (j - 1, i), (j, i + 1), (j, i - 1))
            for jj, ii in neighbours:
                if 0 <= jj < ny and 0 <= ii < nx and mask[jj * nx + ii]:
                    near.append(k)
                    break
    return near


def load_gate():
    """the bracket's own gate, which carries the settlement criterion."""
    gate_path = HERE / "gates.py"
    if not gate_path.exists():
        return None
    spec = importlib.util.spec_from_file_location("wall_bracket_gate", gate_path)
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def leg_case_paths():
    """where each leg's case file lives.

    the cases are kept in the repository, and the run directory holds
    only what the march produced, so the case is sought here rather
    than beside the run's own output.
    """
    return {tag: HERE / f"{tag}.toml" for tag in LEG_TAGS}


LEG_TAGS = (
    "wall-adiabatic",
    "wall-radiating",
    "wall-filmcooled",
    "wall-regenerative",
    "wall-conjugate",
)

BASELINE_TAG = "wall-adiabatic"


def near_wall_series(reader, case, leg_root: Path, window_fields):
    """the near-wall temperature of one leg, field by field, on the window.

    each field contributes its near-wall median and its coldest cell, and
    only the fields in the judged window contribute.
    """
    paths = sorted(glob.glob(str(leg_root / "vtk-*.vtk")))
    out: dict[str, object] = {"field_count": len(paths)}
    if len(paths) < 4:
        out["status"] = "not enough fields"
        return out

    mask = reader.solid_mask(case)
    wanted = set(window_fields) if window_fields else None

    medians = {}
    minima = {}
    for path in paths:
        name = os.path.basename(path)
        digits = ""
        for part in name.split("-"):
            if part.endswith(".vtk"):
                digits = part[:-4]
                break
        try:
            index = int(digits)
        except ValueError:
            continue
        if wanted is not None and index not in wanted:
            continue
        fields, nx, ny, _origin, _spacing = reader.read_vtk(path)
        near = near_wall_cells(mask, nx, ny)
        if not near:
            out["status"] = "the fixture has no gas against the wall"
            return out
        temps = [fields["temperature"][k] for k in near]
        medians[index] = statistics.median(temps)
        minima[index] = min(temps)

    if not medians:
        out["status"] = "the judged window covers none of the saved fields"
        return out

    out["window_fields_used"] = len(medians)
    out["near_wall_mean_k"] = statistics.fmean(list(medians.values()))
    out["near_wall_min_k"] = statistics.fmean(list(minima.values()))
    out["near_wall_coldest_cell_k"] = min(minima.values())
    out["near_wall_mean_spread_k"] = max(medians.values()) - min(medians.values())

    bulk = []
    for index in sorted(medians):
        path = os.path.join(str(leg_root), "vtk-%06d.vtk" % index)
        if not os.path.exists(path):
            continue
        fields, nx, ny, _origin, _spacing = reader.read_vtk(path)
        gas = [fields["temperature"][k] for k in range(nx * ny) if not mask[k]]
        if gas:
            bulk.append(statistics.median(gas))
    if bulk:
        out["bulk_mean_k"] = statistics.fmean(bulk)
        out["near_wall_minus_bulk_k"] = out["near_wall_mean_k"] - out["bulk_mean_k"]

    out["_medians"] = medians
    out["_minima"] = minima
    return out


def paired_departure(baseline, other, key="_medians"):
    """the leg's near-wall departure from the adiabatic leg, field by field."""
    a = baseline.get(key) or {}
    b = other.get(key) or {}
    shared = sorted(set(a) & set(b))
    if not shared:
        return None
    differences = [b[index] - a[index] for index in shared]
    mean = statistics.fmean(differences)
    variance = (
        sum((d - mean) ** 2 for d in differences) / (len(differences) - 1)
        if len(differences) > 1
        else 0.0
    )
    return {
        "delta_k": mean,
        "spread_k": variance ** 0.5,
        "fields": len(shared),
    }


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 2

    run_root = Path(argv[0]).resolve()
    reader = load_field_reader()
    gate = load_gate()
    if gate is None:
        print(json.dumps({"error": "the bracket's gate is missing"}))
        return 1

    # the judged window comes from the gate's own settlement criterion,
    # so every instrument in the bracket speaks about the same stretch of
    # the march.
    rows = []
    for tag in LEG_TAGS:
        leg = gate.measure_one_leg(tag, str(run_root))
        leg["tag"] = tag
        rows.append(leg)
    cross = gate.cross_leg_settlement.assess(rows)
    window_fields = cross.get("window_fields") or []

    cases = leg_case_paths()
    legs = {}
    for tag in LEG_TAGS:
        root = run_root / tag
        case_path = cases[tag]
        if not root.exists() or not case_path.exists():
            continue
        case = reader.read_case(str(case_path))
        legs[tag] = near_wall_series(reader, case, root, window_fields)

    # the adiabatic leg is the baseline: it is the same march without the
    # wall's exchange.
    baseline = legs.get(BASELINE_TAG)
    if isinstance(baseline, dict):
        for tag, entry in legs.items():
            if tag == BASELINE_TAG or not isinstance(entry, dict):
                continue
            if "near_wall_mean_k" not in entry or "near_wall_mean_k" not in baseline:
                continue
            delta = entry["near_wall_mean_k"] - baseline["near_wall_mean_k"]
            entry["near_wall_vs_adiabatic_k"] = delta
            entry["near_wall_colder_than_adiabatic"] = delta < 0.0
            paired = paired_departure(baseline, entry)
            if paired is not None:
                entry["paired_departure"] = paired

    for entry in legs.values():
        if isinstance(entry, dict):
            entry.pop("_medians", None)
            entry.pop("_minima", None)

    print(
        json.dumps(
            {
                "run_root": str(run_root),
                "window_fields": len(window_fields),
                "window_indices": window_fields,
                "settlement_qualified": cross.get("qualified"),
                "settlement_reason": cross.get("reason", ""),
                "legs": legs,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

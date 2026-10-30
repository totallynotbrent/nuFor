#!/usr/bin/env python3
"""the near-wall gas temperature: does the wall reach the gas it touches?

WHY THIS EXISTS

the wall's heat flux has to reach the gas on the other side of it. the
visible symptoms of a wall that does not reach its gas are a believable
reported duty and a gas that does not care, so the interesting
measurement is the gas's response, not the wall's claim.

this module measures the temperature in the cells ADJACENT TO THE WALL
and compares the wall's legs against the adiabatic leg. the near-wall
gas is where the wall's influence must show up first and strongest, and
averaging it over a whole flow period washes out the periodic motion
that hides a weak effect in the bulk averages.

WHAT IT CAUGHT

the viscous axisymmetric lane's radial divergence read the first column's
radial faces for every column. the wall's heat flux was therefore spent
on the gas above the inlet, and the gas under the wall saw nothing. the
reported duty stayed believable while the gas did not respond to it.

USAGE
    python3 near_wall_response.py <run_root> [<leg> ...]
"""

from __future__ import annotations

import glob
import json
import os
import re
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


def leg_profile(reader, case, leg_root: Path) -> dict[str, object]:
    """the near-wall temperature profile of one leg, over its saved fields.

    each field contributes its near-wall mean and minimum, and the leg's
    profile is the median across the fields. the median rather than the
    mean, so one odd field cannot swing the result.
    """
    paths = sorted(glob.glob(str(leg_root / "vtk-*.vtk")))
    out: dict[str, object] = {"field_count": len(paths)}
    if len(paths) < 4:
        out["status"] = "not enough fields"
        return out

    mask = reader.solid_mask(case)
    means = []
    minima = []
    for path in paths:
        fields, nx, ny, origin, spacing = reader.read_vtk(path)
        near = near_wall_cells(mask, nx, ny)
        if not near:
            out["status"] = "the fixture has no gas against the wall"
            return out
        temps = [fields["temperature"][k] for k in near]
        means.append(statistics.median(temps))
        minima.append(min(temps))

    out["near_wall_mean_k"] = statistics.median(means)
    out["near_wall_min_k"] = statistics.median(minima)
    out["near_wall_mean_spread_k"] = (
        statistics.median([abs(m - statistics.median(means)) for m in means])
    )
    out["near_wall_min_spread_k"] = statistics.median(
        [abs(m - statistics.median(minima)) for m in minima]
    )

    # the bulk, for the reference of how much colder the near-wall gas is
    bulk = []
    for path in paths[len(paths) // 2 :]:
        fields, nx, ny, origin, spacing = reader.read_vtk(path)
        gas = [fields["temperature"][k] for k in range(nx * ny) if not mask[k]]
        bulk.append(statistics.median(gas))
    out["bulk_mean_k"] = statistics.median(bulk)
    out["near_wall_minus_bulk_k"] = out["near_wall_mean_k"] - out["bulk_mean_k"]
    return out


def period_averaged_profile(reader, case, leg_root: Path, period: int | None):
    """the near-wall profile of one leg, averaged over its flow period.

    comparing single fields between two marches compares two arbitrary
    phases, and the chamber's periodic motion dominates the difference:
    the raw column deltas swing by hundreds of kelvin in both directions
    while the legs are in fact equivalent. the flow period tiles the
    motion evenly, so the MEAN over one whole period is the cycle's
    average and is free of the phase.

    the fields are selected by their STEP NUMBER, not by their position
    in the file, so two legs that saved different counts still compare
    over the same stretch of the march.

    the period comes from the caller. the bracket's gate detects it, by
    exact field repetition where the orbit has closed and by signal
    recurrence where it has not.
    """
    paths = sorted(glob.glob(str(leg_root / "vtk-*.vtk")))
    out: dict[str, object] = {"field_count": len(paths)}
    if len(paths) < 4:
        out["status"] = "not enough fields"
        return out

    # the step each field was saved at
    stepped = []
    for path in paths:
        digits = re.search(r"(\d+)", os.path.basename(path))
        stepped.append((int(digits.group(1)) if digits else 0, path))
    stepped.sort()

    mask = reader.solid_mask(case)

    def near_mean(step: int, path: str):
        fields, nx, ny, origin, spacing = reader.read_vtk(path)
        near = near_wall_cells(mask, nx, ny)
        if not near:
            return None
        return statistics.mean([fields["temperature"][k] for k in near])

    if period:
        # the last whole period of saved fields
        chosen = stepped[-period:]
    else:
        # the settled half of the march
        chosen = stepped[len(stepped) // 2 :]

    means = []
    for step, path in chosen:
        value = near_mean(step, path)
        if value is not None:
            means.append(value)

    if not means:
        out["status"] = "the fixture has no gas against the wall"
        return out

    out["period"] = period
    out["window_fields"] = len(chosen)
    out["near_wall_mean_k"] = statistics.mean(means)
    out["near_wall_window_spread_k"] = max(means) - min(means)
    return out


def detect_period(run_root: Path):
    """the flow period in saved fields, from the bracket's own gate.

    the gate carries the detection: the checksum criterion where the
    orbit closed exactly, the signal estimate where it did not, and the
    agreement between the legs where neither alone is conclusive.
    reusing it keeps one definition of the period across the bracket's
    instruments.
    """
    import importlib.util

    gate_path = HERE / "gates.py"
    if not gate_path.exists():
        return None, "no gate"
    spec = importlib.util.spec_from_file_location("wall_bracket_gate", gate_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    reader = load_field_reader()
    legs = []
    for tag in ("wall-adiabatic", "wall-radiating", "wall-filmcooled"):
        root = run_root / tag
        case_path = root / "case.toml"
        if not case_path.exists():
            continue
        legs.append(
            {
                "tag": tag,
                "case": reader.read_case(str(case_path)),
                "fields": sorted(glob.glob(str(root / "vtk-*.vtk"))),
            }
        )
    if not legs:
        return None, "no legs"

    period, strength, per_leg, provenance = module.detect_period_across_legs(legs)
    return period, provenance


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 2

    run_root = Path(argv[0]).resolve()
    reader = load_field_reader()

    # the flow period, taken from the bracket's own gate when it has
    # detected one. the gate knows the period by field repetition where
    # the orbit closed and by signal recurrence where it did not.
    period, period_provenance = detect_period(run_root)

    legs = {}
    for tag in ("wall-adiabatic", "wall-radiating", "wall-filmcooled"):
        root = run_root / tag
        if not root.exists():
            continue
        case_path = root / "case.toml"
        if not case_path.exists():
            continue
        case = reader.read_case(str(case_path))
        legs[tag] = period_averaged_profile(reader, case, root, period)

    # the adiabatic leg is the baseline: it is the same march without the
    # wall's exchange.
    baseline = legs.get("wall-adiabatic")
    if isinstance(baseline, dict):
        for tag, entry in legs.items():
            if tag == "wall-adiabatic" or not isinstance(entry, dict):
                continue
            if "near_wall_mean_k" not in entry or "near_wall_mean_k" not in baseline:
                continue
            delta = entry["near_wall_mean_k"] - baseline["near_wall_mean_k"]
            entry["near_wall_vs_adiabatic_k"] = delta
            entry["near_wall_colder_than_adiabatic"] = delta < 0.0

    print(
        json.dumps(
            {
                "run_root": str(run_root),
                "period": period,
                "period_provenance": period_provenance,
                "legs": legs,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""gate finite-rate GRI-Mech source terms on a spatial composition.

run from the repository root:

    python3 cases/spatial-chemistry/reaction_gate.py

the gate runs the public example, which compares one accepted
finite-rate step of the spatial source pass against the reactor's own
implicit step on the same parcel, cell by cell.
"""

import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUTPUT = pathlib.Path(__file__).with_name("reaction_gates.json")


def main():
    """run the gate and write its measurements."""
    run = subprocess.run(
        [
            "cargo",
            "run",
            "-q",
            "-p",
            "nufor-core",
            "--example",
            "spatial_reaction_gate",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    line = [
        candidate
        for candidate in run.stdout.splitlines()
        if candidate.startswith("{")
    ]
    if not line:
        raise SystemExit("the spatial reaction gate produced no json")
    measured = json.loads(line[-1])

    checks = {
        "every_cell_reacted": measured["cells_reacted"] == measured["cells"],
        "nothing_deferred": measured["cells_deferred"] == 0,
        "field_matches_reactor_within_rounding": measured[
            "worst_relative_error_vs_reactor"
        ]
        < 1.0e-9,
        "compositions_stay_normalized": measured["max_sum_error"] < 1.0e-12,
    }
    artifact = {
        "tag": "spatial-reaction",
        "measured": measured,
        "checks": checks,
        "passed": all(checks.values()),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print(json.dumps(artifact, indent=2))
    if not artifact["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""gate per-cell reaction-set selection against the dense reference.

run from the repository root:

    python3 cases/spatial-chemistry/selection_gate.py

the gate runs the public example, which marches a mixed field of
methane-air parcels twice over the same horizon: once with the full
53-species mechanism and once with a per-cell selected subset, and
compares the trajectories.
"""

import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUTPUT = pathlib.Path(__file__).with_name("selection_gates.json")


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
            "spatial_selection_gate",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    line = [c for c in run.stdout.splitlines() if c.startswith("{")]
    if not line:
        raise SystemExit("the selection gate produced no json")
    measured = json.loads(line[-1])

    checks = {
        "selection_actually_shrunk_the_set": measured["mean_active_species"] < 53.0,
        "no_steps_refused_on_either_side": (
            measured["dense_refusals"] == 0 and measured["selected_refusals"] == 0
        ),
        "major_species_track_the_dense_reference": measured[
            "worst_major_relative_divergence"
        ]
        < 1.0e-6,
        "temperature_tracks_the_dense_reference": measured[
            "worst_temperature_relative_divergence"
        ]
        < 1.0e-3,
    }
    artifact = {
        "tag": "spatial-selection",
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

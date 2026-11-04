#!/usr/bin/env python3
"""gate the scheme-width spatial composition seam.

run from the repository root:

    python3 cases/spatial-chemistry/gates.py

the gate builds and runs the public example, then checks that a real
53-specie GRI-Mech composition crosses one planar flow step without
losing a species row or drifting away from a unit composition sum.
"""

import json
import os
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUTPUT = pathlib.Path(__file__).with_name("gates.json")


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
            "spatial_chemistry_gate",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    line = [candidate for candidate in run.stdout.splitlines() if candidate.startswith("{")]
    if not line:
        raise SystemExit("the spatial chemistry gate produced no json")
    measured = json.loads(line[-1])

    checks = {
        "gri30_species_count": measured["species"] == 53,
        "methane_crossed_the_step": measured["transferred_methane_fraction"] > 0.15,
        "all_species_remained_present": measured["minimum_positive_species_fraction"] > 0.0,
        "composition_sum_closed": measured["maximum_composition_sum_error"] < 1.0e-12,
    }
    artifact = {
        "tag": "spatial-chemistry",
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

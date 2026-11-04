#!/usr/bin/env python3
"""gate the blown-in film's injection accounting.

run from the repository root:

    python3 cases/film-injection/gates.py

the gate builds and runs the public example, then checks that the
coolant the chamber received matches the coolant the report published,
under a wall that states no temperature at all.
"""

import json
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
            "film_injection_gate",
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
        raise SystemExit("the film injection gate produced no json")
    measured = json.loads(line[-1])

    checks = {
        "the_adiabatic_wall_admitted_the_coolant": measured["chamber_mass_gain_kg"] > 0.0,
        "the_chamber_took_the_reported_coolant": measured["relative_mass_error"] < 1.0e-12,
        "the_film_met_a_real_area": measured["film_area_m2"] > 0.0,
        "the_film_brought_its_own_composition": measured["received_nitrogen_fraction"] > 0.0,
    }
    artifact = {
        "tag": "film-injection",
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

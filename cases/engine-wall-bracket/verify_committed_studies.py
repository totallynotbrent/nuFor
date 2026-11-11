#!/usr/bin/env python3
"""prove that the committed settlement studies still match the committed record.

the power, search and calibration studies are shipped as evidence: the page
quotes their numbers, and the claim register registers them. they are built
from surrogate draws over the committed record, and every draw is keyed by a
fixed seed, so reproducing them must yield the same bytes.

this script re-runs all three against the record and compares the result with
what is committed. a mismatch means the shipped numbers no longer follow from
the shipped evidence.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile


HERE = os.path.dirname(os.path.abspath(__file__))
RECORD_DIR = os.path.join(HERE, "record")
STUDY_DIR = os.path.join(HERE, "studies")

POWER_PAIRS = 800
SEARCH_PAIRS = 400
CALIBRATION_TRIALS = 400


EXIT_OK = 0
EXIT_MISMATCH = 1
EXIT_ERROR = 2


def canonical(payload):
    """the payload as a comparable blob, independent of how it was indented."""
    return json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()


def run_study(script, argv, root, takes_positional_out=True):
    """run one study against the record and return the payload it wrote."""
    out = os.path.join(root, os.path.basename(script).replace(".py", ".json"))
    environ = dict(os.environ)
    environ.pop("PYTHONPATH", None)
    environ["NFOR_VIEWER"] = "none"
    if takes_positional_out:
        command = [sys.executable, os.path.join(HERE, script)] + argv + [out]
    else:
        command = (
            [sys.executable, os.path.join(HERE, script)] + argv + ["--json", out]
        )
    proc = subprocess.run(
        command, cwd=HERE, env=environ, capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "%s failed with rc=%d: %s" % (script, proc.returncode, proc.stderr[-800:])
        )
    with open(out) as handle:
        return json.load(handle)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--studies",
        default=STUDY_DIR,
        help="the directory holding the committed study artefacts",
    )
    args = parser.parse_args(argv)
    study_dir = os.path.abspath(args.studies)

    if not os.path.isdir(RECORD_DIR):
        print("FAIL: the committed record is missing at %s." % RECORD_DIR)
        return EXIT_ERROR

    expected = {
        "power.json": (
            "detection_power.py",
            [str(POWER_PAIRS), RECORD_DIR],
            True,
        ),
        "search.json": (
            "search_multiplicity.py",
            [str(SEARCH_PAIRS), RECORD_DIR],
            True,
        ),
        "calibration.json": (
            "calibrate_settlement.py",
            [RECORD_DIR, "--trials", str(CALIBRATION_TRIALS)],
            False,
        ),
    }

    root = tempfile.mkdtemp(prefix="verify-studies-")
    mismatches = []
    try:
        for name, (script, argv, positional_out) in sorted(expected.items()):
            committed_path = os.path.join(study_dir, name)
            if not os.path.exists(committed_path):
                print("FAIL: the committed study %s is missing." % committed_path)
                return EXIT_ERROR
            with open(committed_path) as handle:
                committed = json.load(handle)

            produced = run_study(script, argv + ["--record"], root, positional_out)

            if canonical(produced) == canonical(committed):
                print("PASS: %s reproduces from the committed record." % name)
            else:
                mismatches.append(name)
                print("FAIL: %s does not reproduce from the committed record." % name)
    except Exception as exc:
        print("FAIL: the studies could not be re-run: %s" % exc)
        return EXIT_ERROR
    finally:
        shutil.rmtree(root, ignore_errors=True)

    if mismatches:
        print(
            "FAIL: these shipped studies no longer follow from the committed "
            "record: %s" % ", ".join(mismatches)
        )
        return EXIT_MISMATCH
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""prove that the committed gate artefact still matches the committed record."""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile


HERE = os.path.dirname(os.path.abspath(__file__))
CASES_ROOT = os.path.dirname(HERE)
COMMITTED_GATES = os.path.join(HERE, "gates.json")
RECORD_DIR = os.path.join(HERE, "record")


EXIT_OK = 0
EXIT_MISMATCH = 1
EXIT_ERROR = 2
EXIT_NO_ARTEFACT = 3


def digest_bytes(data):
    """the sha256 of some bytes."""
    return hashlib.sha256(data).hexdigest()


def compare(left, right, tolerance=5e-9):
    """where two verdicts differ, treating numbers as equal within tolerance.

    the artefact rounds its numbers before writing them, so a rerun on
    another platform can still differ in the last bits if the platform's
    maths library rounds a borderline case the other way. a number that
    agrees to nine significant decimals is the same number for this
    judgement, so the comparison ignores differences below that and
    reports everything else.
    """
    problems = []

    def walk(one, other, path):
        if isinstance(one, dict) and isinstance(other, dict):
            for key in sorted(set(one) | set(other)):
                if key not in one:
                    problems.append("%s.%s is only in the committed artefact" % (path, key))
                elif key not in other:
                    problems.append("%s.%s is missing from the committed artefact" % (path, key))
                else:
                    walk(one[key], other[key], "%s.%s" % (path, key))
            return
        if isinstance(one, list) and isinstance(other, list):
            if len(one) != len(other):
                problems.append(
                    "%s holds %d entries where the committed artefact holds %d"
                    % (path, len(one), len(other))
                )
                return
            for index, (left_item, right_item) in enumerate(zip(one, other)):
                walk(left_item, right_item, "%s[%d]" % (path, index))
            return
        if isinstance(one, bool) or isinstance(other, bool):
            if one != other:
                problems.append("%s is %r where the committed artefact has %r" % (path, one, other))
            return
        if isinstance(one, (int, float)) and isinstance(other, (int, float)):
            scale = max(abs(float(one)), abs(float(other)), 1.0)
            if abs(float(one) - float(other)) > tolerance * scale:
                problems.append(
                    "%s is %.12g where the committed artefact has %.12g"
                    % (path, float(one), float(other))
                )
            return
        if one != other:
            problems.append(
                "%s is %r where the committed artefact has %r" % (path, one, other)
            )

    walk(left, right, "artefact")
    return problems


def mirror_case(root):
    """build a throwaway tree holding symlinks to the case, its sibling case and the record."""
    sandbox = os.path.join(root, "cases")
    os.makedirs(sandbox)
    case_link = os.path.join(sandbox, "engine-wall-bracket")
    os.mkdir(case_link)
    for name in sorted(os.listdir(HERE)):
        if name in ("gates.json", "__pycache__"):
            continue
        os.symlink(os.path.join(HERE, name), os.path.join(case_link, name))
    sibling = os.path.join(CASES_ROOT, "chamber-throat")
    if os.path.isdir(sibling):
        os.symlink(sibling, os.path.join(sandbox, "chamber-throat"))
    return case_link


def regenerate(root):
    """run the gate against the committed record in a sandbox and return the produced bytes."""
    case_dir = mirror_case(root)
    script = os.path.join(case_dir, "gates.py")
    environ = dict(os.environ)
    environ.pop("PYTHONPATH", None)
    environ["NUFOR_BRACKET_RUNS"] = os.path.join(case_dir, "runs-absent")
    environ["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, script],
        cwd=root,
        env=environ,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    produced = os.path.join(case_dir, "gates.json")
    if not os.path.exists(produced):
        raise RuntimeError(
            "the gate program did not produce an artefact (exit %d)\n%s"
            % (proc.returncode, proc.stdout.strip())
        )
    with open(produced, "rb") as handle:
        return handle.read(), proc.stdout


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gates",
        default=COMMITTED_GATES,
        help="the gate artefact to check against the committed record",
    )
    args = parser.parse_args(argv)
    committed_path = os.path.abspath(args.gates)

    if not os.path.isdir(RECORD_DIR):
        print(
            "FAIL: the committed record is missing at %s, so the gate cannot be re-judged."
            % RECORD_DIR
        )
        return EXIT_ERROR
    manifest = os.path.join(RECORD_DIR, "manifest.json")
    if not os.access(manifest, os.R_OK):
        print(
            "FAIL: the committed record is unreadable at %s, so the gate cannot be re-judged."
            % RECORD_DIR
        )
        return EXIT_ERROR

    if not os.path.exists(committed_path):
        print(
            "MISSING: there is no committed gate artefact yet at %s." % committed_path
        )
        return EXIT_NO_ARTEFACT
    try:
        with open(committed_path, "rb") as handle:
            committed = handle.read()
    except OSError as exc:
        print("FAIL: the committed gate artefact is unreadable: %s" % exc)
        return EXIT_ERROR

    root = tempfile.mkdtemp(prefix="verify-committed-gate-")
    try:
        produced, log = regenerate(root)
    except Exception as exc:
        print("FAIL: the gate program could not re-judge the committed record: %s" % exc)
        return EXIT_ERROR
    finally:
        shutil.rmtree(root, ignore_errors=True)

    if produced == committed:
        print(
            "PASS: the committed engine wall bracket gate artefact matches the "
            "committed record byte-for-byte, sha256=%s (%s)"
            % (digest_bytes(committed), committed_path)
        )
        return EXIT_OK

    # THE BYTES CAN MOVE WHILE THE JUDGEMENT DOES NOT. the maths
    # libraries on different platforms round the last bits of a sum
    # differently, so two runs of the same judgement can disagree in the
    # fifteenth digit. the artefact rounds its numbers before writing
    # them, which makes that rare, but a rerun still deserves a verdict
    # on the numbers themselves rather than a failure on the encoding.
    try:
        with open(committed_path, "r") as handle:
            committed_verdict = json.load(handle)
        produced_verdict = json.loads(produced.decode("utf-8"))
    except (OSError, ValueError) as exc:
        print(
            "FAIL: the committed record and the committed artefact disagree for %s: "
            "regenerated sha256=%s, committed sha256=%s, and the two could not be "
            "compared as verdicts: %s. The shipped evidence is inconsistent, so the "
            "shipped artefact is not the judgement of the committed record."
            % (
                committed_path,
                digest_bytes(produced),
                digest_bytes(committed),
                exc,
            )
        )
        return EXIT_MISMATCH

    differences = compare(produced_verdict, committed_verdict)
    if not differences:
        print(
            "PASS: the committed engine wall bracket gate artefact matches the "
            "committed record, sha256=%s (%s); the two encodings differ but every "
            "quantity agrees."
            % (digest_bytes(committed), committed_path)
        )
        return EXIT_OK

    print(
        "FAIL: the committed record and the committed artefact disagree for %s: "
        "regenerated sha256=%s, committed sha256=%s. The shipped evidence is "
        "inconsistent, so the shipped artefact is not the judgement of the "
        "committed record."
        % (
            committed_path,
            digest_bytes(produced),
            digest_bytes(committed),
        )
    )
    for difference in differences[:20]:
        print("  - %s" % difference)
    if len(differences) > 20:
        print("  - and %d more differences" % (len(differences) - 20))
    return EXIT_MISMATCH


if __name__ == "__main__":
    sys.exit(main())

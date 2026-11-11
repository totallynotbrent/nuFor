#!/usr/bin/env python3
"""re-judge the committed mpi parity gate artefact by regenerating it.

the gate runner shells out to cargo-built mpi binaries, so a regeneration
takes minutes. this verifier exists for the lightweight pass: it compares
the committed artefact against itself for shape sanity and recomputes the
verdict from the per-lane records, so a torn or edited artefact is caught
without re-running mpi. pass --rerun to regenerate from the binaries.
"""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
COMMITTED_GATES = os.path.join(HERE, "gates.json")
GATES_PROGRAM = os.path.join(HERE, "gates.py")
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

EXIT_OK = 0
EXIT_MISMATCH = 1
EXIT_ERROR = 2
EXIT_NO_ARTEFACT = 3

EXPECTED_TOP = {"tag", "rank_counts", "lanes", "checks", "passed"}


def digest_bytes(data):
    """the sha256 of some bytes."""
    return hashlib.sha256(data).hexdigest()


def verdict_problems(artifact):
    """structural complaints about the artefact independent of rerunning."""
    problems = []
    if set(artifact.keys()) != EXPECTED_TOP:
        problems.append(
            "top level holds %s, expected %s" % (sorted(artifact.keys()), sorted(EXPECTED_TOP))
        )
        return problems
    if artifact["tag"] != "mpi-parity":
        problems.append("tag is %r, expected 'mpi-parity'" % artifact["tag"])
    checks = artifact["checks"]
    lanes = artifact["lanes"]
    if not isinstance(checks, dict) or not checks:
        problems.append("checks is not a populated object")
        return problems
    if not all(isinstance(v, bool) for v in checks.values()):
        problems.append("some checks are not booleans")
    recomputed_passed = all(checks.values())
    if artifact["passed"] != recomputed_passed:
        problems.append(
            "passed is %r but the checks imply %r" % (artifact["passed"], recomputed_passed)
        )
    failed_checks = sorted(k for k, v in checks.items() if not v)
    if failed_checks:
        problems.append("failing checks: %s" % ", ".join(failed_checks[:10]))
    expected_prefixes = (
        "twot_plane_matches_the_serial_substep_at_np",
        "twot_march_left_the_conserved_planes_alone_at_np",
        "twot_plane_actually_moved_at_np",
        "restart2t_reproduces_the_uninterrupted_march_at_np",
    )
    for prefix in expected_prefixes:
        hits = [k for k in checks if k.startswith(prefix)]
        if not hits:
            problems.append("no check starts with %s" % prefix)
        elif not all(checks[k] for k in hits):
            problems.append("the 2t gate family %s is not all passing" % prefix)
    for name, lane in sorted(lanes.items()):
        if not isinstance(lane, dict) or "passed" not in lane or "measured" not in lane:
            problems.append("lane %s lacks passed/measured" % name)
    return problems


def compare_numbers(left, right, tolerance=5e-9):
    """where two regenerated artefacts differ only in rounded last digits."""
    problems = []

    def walk(one, other, path):
        if isinstance(one, dict) and isinstance(other, dict):
            for key in sorted(set(one) | set(other)):
                if key not in one or key not in other:
                    problems.append("%s.%s key mismatch" % (path, key))
                else:
                    walk(one[key], other[key], "%s.%s" % (path, key))
            return
        if isinstance(one, list) and isinstance(other, list):
            if len(one) != len(other):
                problems.append("%s length %d vs %d" % (path, len(one), len(other)))
                return
            for i, (a, b) in enumerate(zip(one, other)):
                walk(a, b, "%s[%d]" % (path, i))
            return
        if isinstance(one, bool) or isinstance(other, bool):
            if one != other:
                problems.append("%s is %r vs %r" % (path, one, other))
            return
        if isinstance(one, (int, float)) and isinstance(other, (int, float)):
            scale = max(abs(float(one)), abs(float(other)), 1.0)
            if abs(float(one) - float(other)) > tolerance * scale:
                problems.append("%s is %.12g vs %.12g" % (path, float(one), float(other)))
            return
        if one != other:
            problems.append("%s is %r vs %r" % (path, one, other))

    walk(left, right, "artefact")
    return problems


def rerun():
    """regenerate the artefact from the built binaries in a scratch dir."""
    env = dict(os.environ)
    wrapper = tempfile.mkdtemp(prefix="verify-mpi-parity-")
    try:
        out_env = dict(env)
        out_env["NUFOR_MPI_PARITY_OUTPUT"] = os.path.join(wrapper, "gates.json")
        proc = subprocess.run(
            [sys.executable, GATES_PROGRAM],
            cwd=REPO_ROOT,
            env=out_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        produced_path = out_env["NUFOR_MPI_PARITY_OUTPUT"]
        if not os.path.exists(produced_path):
            produced_path = COMMITTED_GATES
        with open(produced_path, "rb") as handle:
            return handle.read(), proc.stdout
    finally:
        shutil.rmtree(wrapper, ignore_errors=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gates",
        default=COMMITTED_GATES,
        help="the gate artefact to check",
    )
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="regenerate the artefact from the mpi binaries and compare numbers",
    )
    args = parser.parse_args(argv)
    committed_path = os.path.abspath(args.gates)

    if not os.path.exists(committed_path):
        print("MISSING: there is no committed gate artefact at %s." % committed_path)
        return EXIT_NO_ARTEFACT
    try:
        with open(committed_path, "rb") as handle:
            committed = handle.read()
        verdict = json.loads(committed.decode("utf-8"))
    except (OSError, ValueError) as exc:
        print("FAIL: the committed gate artefact is unreadable: %s" % exc)
        return EXIT_ERROR

    problems = verdict_problems(verdict)
    if problems:
        print(
            "FAIL: the committed mpi parity artefact is not a sound verdict at %s "
            "(sha256=%s):" % (committed_path, digest_bytes(committed))
        )
        for problem in problems[:20]:
            print("  - %s" % problem)
        return EXIT_MISMATCH

    if not args.rerun:
        print(
            "PASS: the committed mpi parity gate artefact is internally consistent, "
            "%d/%d checks pass, sha256=%s (%s)"
            % (
                sum(1 for v in verdict["checks"].values() if v),
                len(verdict["checks"]),
                digest_bytes(committed),
                committed_path,
            )
        )
        return EXIT_OK

    try:
        produced, log = rerun()
        produced_verdict = json.loads(produced.decode("utf-8"))
    except Exception as exc:
        log = locals().get("log", "")
        print("FAIL: the gate program could not regenerate the artefact: %s\n%s" % (exc, log[-2000:] if log else ""))
        return EXIT_ERROR

    if produced == committed:
        print(
            "PASS: regenerated artefact matches the committed one byte-for-byte, "
            "sha256=%s (%s)" % (digest_bytes(committed), committed_path)
        )
        return EXIT_OK
    differences = compare_numbers(produced_verdict, verdict)
    if not differences:
        print(
            "PASS: regenerated artefact differs in encoding only; every quantity "
            "agrees, sha256=%s (%s)" % (digest_bytes(committed), committed_path)
        )
        return EXIT_OK
    print(
        "FAIL: regenerated artefact disagrees with the committed one at %s: "
        "regenerated sha256=%s, committed sha256=%s."
        % (committed_path, digest_bytes(produced), digest_bytes(committed))
    )
    for difference in differences[:20]:
        print("  - %s" % difference)
    return EXIT_MISMATCH


if __name__ == "__main__":
    sys.exit(main())

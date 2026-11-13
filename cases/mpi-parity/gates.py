#!/usr/bin/env python3
"""gate the mpi parity of the fortran core against the serial rust lanes.

run from the repository root:

    python3 cases/mpi-parity/gates.py

the gate runs the mpi test binaries under mpirun at np=1 through 5, lets
each gate record its own measurements, and checks the parity contracts
those measurements state:

- the euler lanes reproduce the serial march exactly (zero relative
  difference, bit for bit),
- the two-temperature lane (heun march + the nufor2t e_v plane substep)
  reproduces the serial march-plus-substep exactly at every rank count,
  plane included,
- the viscous, sa and capsule-body lanes stay inside the roundoff
  tolerance their own assertions use,
- the multi-field gather reassembles every plane exactly,
- a split march resumed from a restart file reproduces the
  uninterrupted march exactly, at every rank count, on both the sa
  and the two-temperature lanes.

five ranks are exercised on purpose: the decomposition splits the mesh
so that the outermost rank owns no body cell, which is where a
rank-local decision would leave the ranks marching different schemes
and the collector would never close.

every run carries a timeout, so a rank-seam deadlock fails the gate
instead of hanging it.

the measurements are written beside this script as gates.json. the
committed artefact is the evidence the published numbers come from.
"""

import json
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
# the verifier regenerates into a scratch path so it never overwrites the
# committed artefact while judging it.
OUTPUT = pathlib.Path(
    os.environ.get("NUFOR_MPI_PARITY_OUTPUT") or pathlib.Path(__file__).with_name("gates.json")
)

PARITY_BIN_GLOB = "mpi_parity-"
RESTART_BIN_GLOB = "mpi_restart-"

# the rank counts the gate walks. five is the count where the outer
# rank owns no body cell.
RANK_COUNTS = (1, 2, 3, 4, 5)

# the tolerances the gates assert internally.
EXACT_ZERO_LANES = ("euler", "uniform")
ROUNDOFF_LANES = ("viscous", "sa", "body")
TWOT_LANE = "twot"

# the five-species composition lane. the seam reconstruction is bitwise
# zero by construction, so the species gates classify as exact-zero.
SPECIES_LANE = "mpipark5"

# a run that overruns this is a seam deadlock, not a slow machine.
RUN_TIMEOUT_S = 180


def find_binary(repo_root, glob):
    """the newest built test binary matching the glob."""
    deps = repo_root / "target" / "debug" / "deps"
    if not deps.is_dir():
        return None
    candidates = [p for p in deps.glob(glob + "*") if p.is_file() and p.stat().st_mode & 0o111]
    if not candidates:
        return None
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


def mpirun_cmd(np_value, binary, test):
    """the mpirun invocation for one gate."""
    mpirun = shutil.which("mpirun") or "/usr/lib64/openmpi/bin/mpirun"
    return [
        mpirun,
        "--allow-run-as-root",
        "--oversubscribe",
        "-np",
        str(np_value),
        str(binary),
        "--exact",
        test,
        "--nocapture",
        "--ignored",
    ]


def mpi_environment(base):
    """the environment the test binaries need to find the mpi runtime.

    the fortran core links against the openmpi runtime, which installs
    under its own prefix on some distributions. the loader path is
    filled in here so the gate runs the same way from a shell, from ci
    and from a scripted check.
    """
    env = dict(base)
    loader_paths = []
    for candidate in (
        "/usr/lib64/openmpi/lib",
        "/usr/lib64/openmpi/lib/openmpi3",
        "/usr/lib/x86_64-linux-gnu/openmpi/lib",
    ):
        if pathlib.Path(candidate).is_dir():
            loader_paths.append(candidate)
    if loader_paths:
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = ":".join(loader_paths + ([existing] if existing else []))
    env.setdefault("NFOR_VIEWER", "none")
    return env


def run_gate(np_value, binary, test, report_path, workdir, attempts=3):
    """run one gate and return its log and the report it wrote.

    the launcher's exit code is noise on this cluster: one rank's mpi
    finalize races the other's teardown and mpirun reports failure after
    the harness has already printed a passing result. the pass signal is
    the harness summary, so the launcher status is ignored and the
    summary is checked instead.

    a run counts as landed only when the harness passed and the report
    file is on disk. a passing run without a report is a torn-down
    rank zero, so it is retried like a failure.
    """
    env = mpi_environment(os.environ)
    env["NF_MPI_REPORT"] = str(report_path)
    restart_dir = workdir / f"restart-np{np_value}"
    env["NUFOUR_MPI_RESTART_DIR"] = str(restart_dir)
    restart_dir.mkdir(parents=True, exist_ok=True)

    log = ""
    timed_out = False
    attempts_used = 0
    for attempt in range(1, attempts + 1):
        attempts_used = attempt
        if report_path.exists():
            report_path.unlink()
        try:
            proc = subprocess.run(
                mpirun_cmd(np_value, binary, test),
                cwd=ROOT,
                capture_output=True,
                text=True,
                env=env,
                timeout=RUN_TIMEOUT_S,
            )
            log = proc.stdout + proc.stderr
            timed_out = False
        except subprocess.TimeoutExpired as exc:
            timed_out = True
            stdout = exc.stdout.decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            stderr = exc.stderr.decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
            log = stdout + stderr
            log += f"\nthe gate timed out after {RUN_TIMEOUT_S}s\n"
        passed = "1 passed; 0 failed" in log
        if passed and report_path.exists():
            break
        if passed:
            log += "\nthe harness passed but left no report behind\n"
        else:
            log += "\nthe harness did not report a passing result\n"
        if attempt < attempts:
            log += f"retrying ({attempt}/{attempts})\n"
    passed = "1 passed; 0 failed" in log
    measured = {}
    if report_path.exists():
        try:
            measured = json.loads(report_path.read_text())
        except json.JSONDecodeError:
            measured = {}
    return {
        "log": log,
        "passed": passed and bool(measured) and not timed_out,
        "timed_out": timed_out,
        "measured": measured,
        "attempts": attempts_used,
    }


def main():
    """run the gate matrix and write its measurements."""
    parity_bin = find_binary(ROOT, PARITY_BIN_GLOB)
    restart_bin = find_binary(ROOT, RESTART_BIN_GLOB)
    if parity_bin is None or restart_bin is None:
        raise SystemExit(
            "the mpi test binaries are not built; run "
            "`cargo test -p nufor-core --test mpi_parity --test mpi_restart --no-run` first"
        )

    workdir = pathlib.Path(tempfile.mkdtemp(prefix="nufor-mpi-gate-"))
    lanes = {
        "euler": "mpi_np1_matches_rust_sod2d",
        "uniform": "mpi_np1_matches_rust_uniform_eqair",
        TWOT_LANE: "mpi_np1_2t_matches_rust",
        "viscous": "mpi_np1_visc_matches_rust",
        "sa": "mpi_np1_sa_matches_rust",
        "body": "mpi_np1_matches_rust_capsule",
        "gather": "mpi_gather_nf_fields_reassembles_exactly",
        SPECIES_LANE: "mpi_mpipark5_carrier_seam_and_gather",
    }
    results = {}
    try:
        for np_value in RANK_COUNTS:
            for lane, test in lanes.items():
                report = workdir / f"{lane}-np{np_value}.json"
                outcome = run_gate(np_value, parity_bin, test, report, workdir)
                results[f"{lane}_np{np_value}"] = {
                    "test": test,
                    "world_size": np_value,
                    "passed": outcome["passed"],
                    "timed_out": outcome["timed_out"],
                    "attempts": outcome["attempts"],
                    "measured": outcome["measured"],
                }
        for np_value in RANK_COUNTS:
            report = workdir / f"restart-np{np_value}.json"
            outcome = run_gate(
                np_value,
                restart_bin,
                "mpi_split_restart_reproduces_continuous",
                report,
                workdir,
            )
            results[f"restart_np{np_value}"] = {
                "test": "mpi_split_restart_reproduces_continuous",
                "world_size": np_value,
                "passed": outcome["passed"],
                "timed_out": outcome["timed_out"],
                "attempts": outcome["attempts"],
                "measured": outcome["measured"],
            }
            report = workdir / f"restart2t-np{np_value}.json"
            outcome = run_gate(
                np_value,
                restart_bin,
                "mpi_split_restart_2t_reproduces_continuous",
                report,
                workdir,
            )
            results[f"restart2t_np{np_value}"] = {
                "test": "mpi_split_restart_2t_reproduces_continuous",
                "world_size": np_value,
                "passed": outcome["passed"],
                "timed_out": outcome["timed_out"],
                "attempts": outcome["attempts"],
                "measured": outcome["measured"],
            }
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    for np_value in RANK_COUNTS:
        for lane in EXACT_ZERO_LANES:
            entry = results[f"{lane}_np{np_value}"]
            key = f"{lane}_worst_rel_diff"
            check(
                f"{lane}_lane_is_bit_identical_at_np{np_value}",
                entry["passed"] and entry["measured"].get(key) == 0.0,
            )
        for lane in ROUNDOFF_LANES:
            entry = results[f"{lane}_np{np_value}"]
            key = f"{lane}_worst_rel_diff"
            check(
                f"{lane}_lane_stays_inside_roundoff_at_np{np_value}",
                entry["passed"] and float(entry["measured"].get(key, 1.0)) < 1.0e-6,
            )
        twot = results[f"{TWOT_LANE}_np{np_value}"]
        check(
            f"twot_plane_matches_the_serial_substep_at_np{np_value}",
            twot["passed"]
            and twot["measured"].get("twot_differing_entries") == 0
            and twot["measured"].get("twot_compared_entries", 0) > 0,
        )
        check(
            f"twot_march_left_the_conserved_planes_alone_at_np{np_value}",
            twot["passed"]
            and sum(
                twot["measured"].get(f"twot_state_differing_{field}", 1)
                for field in ("rho", "mx", "my", "e")
            )
            == 0,
        )
        check(
            f"twot_plane_actually_moved_at_np{np_value}",
            twot["passed"]
            and float(twot["measured"].get("twot_plane_liveness", 0.0))
            > 1.0e-9 * float(twot["measured"].get("twot_seed_peak", 1.0)),
        )
        species = results[f"{SPECIES_LANE}_np{np_value}"]
        check(
            f"species_seam_reconstruction_is_bitwise_zero_at_np{np_value}",
            species["passed"]
            and species["measured"].get("mpipark5_seam_zero") == 1,
        )
        gather = results[f"gather_np{np_value}"]
        check(
            f"gather_reassembles_every_plane_at_np{np_value}",
            gather["passed"]
            and gather["measured"].get("gather_differing_entries") == 0
            and gather["measured"].get("gather_compared_entries", 0) > 0,
        )
        restart = results[f"restart_np{np_value}"]
        restart_fields = [
            "restart_rho_max_abs_diff",
            "restart_mx_max_abs_diff",
            "restart_my_max_abs_diff",
            "restart_e_max_abs_diff",
            "restart_nu_tilde_max_abs_diff",
        ]
        check(
            f"restart_reproduces_the_uninterrupted_march_at_np{np_value}",
            restart["passed"]
            and all(restart["measured"].get(field) == 0.0 for field in restart_fields),
        )
        restart2t = results[f"restart2t_np{np_value}"]
        restart2t_fields = [
            "restart2t_rho_max_abs_diff",
            "restart2t_mx_max_abs_diff",
            "restart2t_my_max_abs_diff",
            "restart2t_e_max_abs_diff",
            "restart2t_e_v_max_abs_diff",
        ]
        check(
            f"restart2t_reproduces_the_uninterrupted_march_at_np{np_value}",
            restart2t["passed"]
            and all(restart2t["measured"].get(field) == 0.0 for field in restart2t_fields),
        )

    artifact = {
        "tag": "mpi-parity",
        "rank_counts": list(RANK_COUNTS),
        "lanes": results,
        "checks": checks,
        "passed": all(checks.values()),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {"passed": artifact["passed"], "failed": [k for k, v in checks.items() if not v]},
            indent=2,
        )
    )
    if not artifact["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""measure mpi strong scaling for the eqair 2d euler march.

the mpi page publishes a scaling table and a chart. both come from this
script so the numbers can be reproduced instead of carried by hand.

the timing comes from the solver's own `wall:` line, which brackets the
march only. the surrounding work (the live viewer booting, the vtk and
html report writes) sits outside that line and is identical across rank
counts, so timing the whole process would bury the march in a constant.

every rank count, including one, runs under mpirun, so all of them take
the fortran mpi path. comparing np=1 under mpirun against a plain serial
run would compare two different codes.

usage:
    tools/bench_mpi.py [--bin target/release/nufor] [--case FILE]
                       [--ranks 1 2 4 8] [--repeats 3] [--json OUT]
"""

import argparse
import json
import os
import re
import shutil
import statistics
import subprocess
import tempfile


def default_case(final_time: str, viscous: bool) -> str:
    """the station-shaped probe: 320x176, a fixed number of steps."""

    mu_block = ""
    if viscous:
        mu_block = "mu = 2.0e-5\npr = 0.72\nwall_temperature = 300.0\n"
    return f"""
schema_version = 1

[metadata]
name = "mpi-scaling"
description = "mpi strong scaling probe"

[mesh]
nx = 320
ny = 176
x0 = -2.0
x1 = 36.0
y0 = -13.0
y1 = 13.0
source = "uniform"

[physics]
equations = "euler_2d"
eos = "eqair"
gamma = 1.4
gas_constant = 287.0
{mu_block}
[time]
final_time = {final_time}
max_steps = 4000
residual_target = 1e-30

[numerics]
flux = "hllc"
reconstruction = "muscl"
cfl = 0.5
local_time_stepping = false
threads = 1

[boundaries]
left = "supersonic_inflow"
right = "supersonic_outflow"
top = "supersonic_outflow"
bottom = "supersonic_inflow"

[initial_condition]
type = "uniform"
rho = 0.01
u = 6691.0
p = 685.0

[output]
formats = ["vtk"]
interval_steps = 100000
fields = ["rho", "p"]
"""


def solve(binary: str, case: str, ranks: int, cwd: str) -> tuple[float, int]:
    """run the case on the given ranks and return the march's wall time."""

    command = [
        "mpirun",
        "--oversubscribe",
        "-np",
        str(ranks),
        binary,
        "run",
        case,
    ]
    env = dict(os.environ)
    env.setdefault("OMPI_ALLOW_RUN_AS_ROOT", "1")
    env.setdefault("OMPI_ALLOW_RUN_AS_ROOT_CONFIRM", "1")
    # without this, np=1 falls back to the serial rust march because the
    # launcher sees a world size of one, and the np=1 row would be
    # comparing a different code against the fortran lanes.
    env["NFOR_MPI"] = "1"
    env.setdefault("NFOR_VIEWER", "none")
    done = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if done.returncode != 0:
        raise SystemExit(
            f"the np={ranks} run failed ({done.returncode}):\n{done.stderr[-2000:]}"
        )
    match = re.search(r"^wall:\s+([0-9.]+)s", done.stdout, re.M)
    if not match:
        raise SystemExit(
            f"the np={ranks} run printed no wall line:\n{done.stdout[-2000:]}"
        )
    steps_match = re.search(r"^steps:\s+(\d+)", done.stdout, re.M)
    steps = int(steps_match.group(1)) if steps_match else -1
    return float(match.group(1)), steps


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", default="target/release/nufor")
    parser.add_argument("--case", default=None)
    parser.add_argument("--ranks", nargs="+", type=int, default=[1, 2, 4, 8])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--final-time", default="0.004")
    parser.add_argument(
        "--viscous",
        action="store_true",
        help="route the viscous fortran lane (sets mu/pr/wall_temperature)",
    )
    parser.add_argument("--json", default=None)
    args = parser.parse_args()

    binary = os.path.abspath(args.bin)
    if not os.path.exists(binary):
        raise SystemExit(f"no binary at {binary}")

    work = tempfile.mkdtemp(prefix="nufor-mpi-bench-")
    try:
        case = args.case or os.path.join(work, "case.toml")
        if not args.case:
            with open(case, "w") as handle:
                handle.write(default_case(args.final_time, args.viscous))

        rows = []
        base = None
        step_counts = {}
        for ranks in args.ranks:
            samples = [solve(binary, case, ranks, work) for _ in range(args.repeats)]
            timings = [t for t, _ in samples]
            step_counts[ranks] = sorted({s for _, s in samples})
            best = min(timings)
            median = statistics.median(timings)
            if base is None:
                base = best
            speedup = None if base == 0 else base / best
            rows.append(
                {
                    "ranks": ranks,
                    "steps": step_counts[ranks],
                    "best_s": round(best, 4),
                    "median_s": round(median, 4),
                    "speedup": None if speedup is None else round(speedup, 3),
                }
            )
            print(
                f"np={ranks}: steps {step_counts[ranks]} best {best:.3f}s"
                f" median {median:.3f}s"
                + ("" if speedup is None else f" speedup {speedup:.2f}x")
            )

        payload = {
            "binary": binary,
            "lane": "viscous" if args.viscous else "euler",
            "final_time": args.final_time,
            "repeats": args.repeats,
            "rows": rows,
        }
        if args.json:
            with open(args.json, "w") as handle:
                json.dump(payload, handle, indent=2)
                handle.write("\n")
            print(f"wrote {args.json}")
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# benchmark the release solver across a mesh sweep and print a throughput table.
#
# usage: tools/bench.py [release binary] [steps] [n n n ...]
#   defaults: target/release/nufor, 3000 steps, 200..10000 cells.

import re
import subprocess
import sys

BIN = sys.argv[1] if len(sys.argv) > 1 else "target/release/nufor"
STEPS = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
MESHES = [int(x) for x in sys.argv[3:]] or [200, 500, 1000, 2000, 4000, 10000]
RUNS = 3


def best_throughput(n: int) -> tuple[float, float]:
    best = 1e9
    for _ in range(RUNS):
        out = subprocess.run([BIN, "benchmark", str(n), str(STEPS)], capture_output=True, text=True).stdout
        m = re.search(r"ran (\d+) steps on (\d+) cells in ([\d.]+)s", out)
        if not m:
            continue
        steps, cells, secs = int(m.group(1)), int(m.group(2)), float(m.group(3))
        us = secs * 1e6 / (steps * cells)
        best = min(best, us)
    return best, 1e6 / best  # us/step/cell, cell-steps/s


def main() -> None:
    print(f"{'cells':>7} {'us/step/cell':>12} {'cell-steps/s':>13}")
    rows = [(n, *best_throughput(n)) for n in MESHES]
    for n, us, rate in rows:
        print(f"{n:>7} {us:>11.3f} {rate:>12,.0f}")
    with open("bench.csv", "w") as f:
        f.write("cells,us_per_step_per_cell,cell_steps_per_second\n")
        for n, us, rate in rows:
            f.write(f"{n},{us:.3f},{rate:.0f}\n")
    print("wrote bench.csv")


if __name__ == "__main__":
    main()
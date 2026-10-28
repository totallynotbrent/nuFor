# the catalytic-wall accommodation sweep: the wall heat flux the march
# reports along the surface for each accommodation coefficient, distilled
# from the run logs into samples.csv and summary.csv, so the published
# series can be audited from a clean checkout.
#
#   python3 samples.py            # writes samples.csv and summary.csv
#
# the run logs are the committed record of the runs; the vtk files and
# the case files stay out of the repository.
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))

ACCOMMODATION = ("0", "0.25", "0.5", "0.75", "1")


def read_log(label):
    """one run log: the peak, the continuum probe and the surface points."""
    path = os.path.join(HERE, f"run-log-eta-{label}.txt")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}")
    text = open(path).read()

    def grab(pattern, cast=float):
        m = re.search(pattern, text)
        return cast(m.group(1)) if m else None

    peak = grab(r"peak wall heat flux ([-\d.eE+]+) MW/m\^2")
    peak_x = grab(r"peak wall heat flux [-\d.eE+]+ MW/m\^2 at surface point \(([-\d.eE+]+)")
    peak_y = grab(r"peak wall heat flux [-\d.eE+]+ MW/m\^2 at surface point \([-\d.eE+]+, ([-\d.eE+]+)")
    steps = grab(r"time [-\d.eE+]+ s in (\d+) steps", int)
    wall_s = grab(r"time [-\d.eE+]+ s in \d+ steps, wall ([-\d.eE+]+) s")
    t_max = grab(r"max temperature ([-\d.eE+]+) K")
    kn_rho = grab(r"max kn_gll\(rho\) ([-\d.eE+]+)")
    kn_t = grab(r"max kn_gll\(t\) ([-\d.eE+]+)")
    above = grab(r"cells above breakdown threshold ([-\d.eE+]+) pct")
    m = re.search(
        r"wall heat flux by surface point \(MW/m\^2\):\n((?:\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\n)+)",
        text,
    )
    points = []
    if m:
        for line in m.group(1).strip().splitlines():
            cols = line.split()
            points.append((float(cols[0]), float(cols[1]), float(cols[2])))
    return {
        "eta": float(label),
        "peak_q_mwm2": peak,
        "peak_x": peak_x,
        "peak_y": peak_y,
        "steps": steps,
        "wall_seconds": wall_s,
        "max_temperature_k": t_max,
        "kn_gll_rho": kn_rho,
        "kn_gll_t": kn_t,
        "cells_above_breakdown_pct": above,
        "points": points,
    }


def main():
    runs = [read_log(label) for label in ACCOMMODATION]
    with open(os.path.join(HERE, "samples.csv"), "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["eta", "x", "y", "q_w_mwm2"])
        for run in runs:
            for x, y, q in run["points"]:
                writer.writerow([f"{run['eta']:g}", f"{x:.4f}", f"{y:.4f}", f"{q:.6f}"])
    with open(os.path.join(HERE, "summary.csv"), "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow([
            "eta",
            "peak_q_mwm2",
            "peak_x",
            "peak_y",
            "steps",
            "wall_seconds",
            "max_temperature_k",
            "kn_gll_rho",
            "kn_gll_t",
            "cells_above_breakdown_pct",
        ])
        for run in runs:
            writer.writerow([
                f"{run['eta']:g}",
                f"{run['peak_q_mwm2']:.6f}",
                f"{run['peak_x']:.4f}",
                f"{run['peak_y']:.4f}",
                run["steps"],
                f"{run['wall_seconds']:.1f}",
                f"{run['max_temperature_k']:.1f}",
                f"{run['kn_gll_rho']:.4g}",
                f"{run['kn_gll_t']:.4g}",
                f"{run['cells_above_breakdown_pct']:.2f}",
            ])
    base = runs[0]["peak_q_mwm2"]
    for run in runs:
        amp = run["peak_q_mwm2"] / base if base else float("nan")
        print(f"eta {run['eta']:g}: peak {run['peak_q_mwm2']:.4f} MW/m2 "
              f"(x {run['peak_x']:.3f}, y {run['peak_y']:.3f}), "
              f"amplification {amp:.4f}")


if __name__ == "__main__":
    main()

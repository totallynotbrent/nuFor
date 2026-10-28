# the catalytic-wall accommodation figures: the wall heat flux along
# the surface for each accommodation coefficient, and the peak against
# the accommodation itself. reads the committed sweep samples.
#
#   python3 accommodation_figures.py
#
# the samples come from samples.csv, which samples.py distils from the
# committed run logs; nothing here re-reads a run product.
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def read_samples(path):
    points = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            eta = float(row["eta"])
            points.setdefault(eta, []).append(
                (float(row["x"]), float(row["y"]), float(row["q_w_mwm2"]))
            )
    return points


def sweep_angle(xs, ys):
    # the flow arrives from the left, so the stagnation point is at the
    # upstream end of the upper semicircle; the sweep angle runs from
    # there over the shoulder to the trailing edge
    ang = np.degrees(np.arctan2(ys, xs))
    return 180.0 - ang


def main():
    points = read_samples(os.path.join(HERE, "samples.csv"))
    etas = sorted(points)
    colours = [GRAY, ORANGE, PURPLE, GREEN, RED]

    fig, ax = figure(figsize=(7.5, 3.6))
    for colour, eta in zip(colours, etas):
        pts = np.array(points[eta])
        ang = sweep_angle(pts[:, 0], pts[:, 1])
        ax.plot(ang, pts[:, 2], color=colour, lw=1.8, marker="o", ms=3.0,
                label=f"accommodation {eta:g}")
    ax.set_xlim(0, 180)
    ax.set_xticks(list(range(0, 181, 30)))
    ax.set_xlabel("surface sweep angle from stagnation (deg)")
    ax.set_ylabel(r"$q_w$ (MW/m$^2$)")
    ax.set_title("wall heat flux against the wall accommodation, Mach 20 N2 cylinder")
    ax.legend(fontsize=9)
    save(fig, "cyl-park5visc-accommodation-qw", out_dir=HERE)

    peaks = []
    for eta in etas:
        pts = np.array(points[eta])
        peaks.append((eta, float(pts[:, 2].max())))
    fig2, ax2 = figure(figsize=(7, 3.6))
    ax2.plot([p[0] for p in peaks], [p[1] for p in peaks], color=ORANGE, lw=2.0,
             marker="o", ms=5)
    for eta, peak in peaks:
        ax2.annotate(f"{peak:.2f}", (eta, peak), textcoords="offset points",
                     xytext=(0, 8), ha="center", fontsize=9)
    ax2.set_xlabel("wall accommodation coefficient")
    ax2.set_ylabel(r"peak $q_w$ (MW/m$^2$)")
    lo = min(p[1] for p in peaks)
    hi = max(p[1] for p in peaks)
    ax2.set_ylim(lo - 0.05 * (hi - lo), hi + 0.12 * (hi - lo))
    ax2.set_title("the peak wall heat flux against the wall accommodation")
    save(fig2, "cyl-park5visc-accommodation-peak", out_dir=HERE)


if __name__ == "__main__":
    main()

# park5 viscous catalytic-wall figures on the shared house template:
# the wall heat flux along the surface for the non-catalytic vs fully
# catalytic pair, and the near-wall composition the catalytic condition
# drains. regenerates from the two runs' vtk and their run logs.
import sys
import os
import re
import math

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, GRAY

import numpy as np


def read_fields(path):
    txt = open(path).read()
    blocks = re.findall(
        r"SCALARS (\S+) double 1\nLOOKUP_TABLE default\n((?:[-\d.eE+ ]+\n)+)",
        txt,
    )
    return {name: np.array([float(x) for x in body.split()]) for name, body in blocks}


def read_wall_flux(path):
    # the run log's surface-point block: triples of (x, y, q MW/m2)
    txt = open(path).read()
    m = re.search(
        r"wall heat flux by surface point \(MW/m\^2\):\n((?:\s*[-\d.]+\s+[-\d.]+\s+[-\d.]+\n)+)",
        txt,
    )
    rows = [line.split() for line in m.group(1).strip().splitlines()]
    pts = np.array([[float(c) for c in r[:3]] for r in rows])
    return pts[:, 0], pts[:, 1], pts[:, 2]


def sweep_angle(xs, ys):
    # the flow arrives from the left, so the stagnation point is at the
    # upstream end of the upper semicircle; the sweep angle runs from
    # there over the shoulder to the trailing edge
    ang = np.degrees(np.arctan2(ys, xs))
    return 180.0 - ang


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    f_none = read_fields(os.path.join(here, "cylinder-m20-park5-visc-none.vtk"))
    f_full = read_fields(os.path.join(here, "cylinder-m20-park5-visc-full.vtk"))
    x_none, y_none, q_none = read_wall_flux(os.path.join(here, "run-log-none.txt"))
    x_full, y_full, q_full = read_wall_flux(os.path.join(here, "run-log-full.txt"))

    sw_none = sweep_angle(x_none, y_none)
    sw_full = sweep_angle(x_full, y_full)

    fig, ax = figure(figsize=(7.5, 3.6))
    ax.plot(sw_none, q_none, color=GRAY, lw=1.8, marker="o", ms=3.5, label="non-catalytic")
    ax.plot(sw_full, q_full, color=ORANGE, lw=1.8, marker="o", ms=3.5, label="fully-catalytic")
    ax.set_xlim(0, 180)
    ax.set_xticks(list(range(0, 181, 30)))
    ax.set_xlabel("surface sweep angle from stagnation (deg)")
    ax.set_ylabel(r"$q_w$ (MW/m$^2$)")
    ax.set_title("wall heat flux along the surface, Mach 20 N2 cylinder, 1000 K wall")
    ax.legend()
    save(fig, "cyl-park5visc-qw-bracket", out_dir=here)

    # near-wall atom content: the first fluid row above the symmetry
    # axis is the stagnation row; compare the y_N profile along the
    # stagnation line between the two conditions
    nx, ny = 400, 200
    x0, x1, y0, y1 = -1.8, 5.0, 0.0, 3.4
    dx = (x1 - x0) / nx

    def stagline(f, name):
        xs, ys = [], []
        j = 0
        for i in range(nx):
            x = x0 + (i + 0.5) * dx
            if -1.6 < x < -0.7:
                xs.append(x)
                ys.append(f[name][i + j * nx])
        return np.array(xs), np.array(ys)

    xs, yn_none = stagline(f_none, "y_N")
    _, yn_full = stagline(f_full, "y_N")
    fig2, ax2 = figure(figsize=(7, 3.6))
    ax2.plot(xs, yn_none, color=GRAY, lw=2.0, label="non-catalytic")
    ax2.plot(xs, yn_full, color=PURPLE, lw=2.0, label="fully-catalytic")
    ax2.set_xlabel("x (m)")
    ax2.set_ylabel(r"$Y_\mathrm{N}$ along the stagnation line")
    ax2.set_title("the catalytic wall drains atoms from the near-wall layer")
    ax2.legend()
    save(fig2, "cyl-park5visc-yn-stagline", out_dir=here)


if __name__ == "__main__":
    main()

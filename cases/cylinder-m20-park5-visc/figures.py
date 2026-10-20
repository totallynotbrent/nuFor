# park5 viscous catalytic-wall bracket figures on the shared house
# template: the wall heat flux along the surface for the non-catalytic
# vs fully-catalytic pair, and the near-wall composition the catalytic
# condition drains. regenerates from the two runs' vtk.
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
        r'SCALARS (\S+) double 1\nLOOKUP_TABLE default\n((?:[-\d.eE+ ]+\n)+)',
        txt,
    )
    return {name: np.array([float(x) for x in body.split()]) for name, body in blocks}


def read_wall_flux(path):
    # the run-log's "wall heat flux vs surface x" block: pairs of
    # (surface x, q MW/m2)
    txt = open(path).read()
    m = re.search(r'wall heat flux vs surface x \(MW/m\^2\):\n((?:\s*[-\d.]+\s+[-\d.]+\n)+)', txt)
    pairs = [line.split() for line in m.group(1).strip().splitlines()]
    return np.array([float(p[0]) for p in pairs]), np.array([float(p[1]) for p in pairs])


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    f_none = read_fields(os.path.join(here, "cylinder-m20-park5-visc-none.vtk"))
    f_full = read_fields(os.path.join(here, "cylinder-m20-park5-visc-full.vtk"))
    x_none, q_none = read_wall_flux(os.path.join(here, "run-log-none.txt"))
    x_full, q_full = read_wall_flux(os.path.join(here, "run-log-full.txt"))

    # surface angle from the surface x: the surface points parametrize
    # the upper half of the r=1 circle; theta = acos(x) for the front
    # face, wrapping around the back
    def theta_of(xs):
        return np.array([math.acos(max(-1.0, min(1.0, x))) for x in xs])

    th_none, th_full = theta_of(x_none), theta_of(x_full)

    fig, ax = figure(figsize=(7.5, 3.6))
    ax.plot(th_none * 180 / math.pi, q_none, color=GRAY, lw=2.0, label="non-catalytic")
    ax.plot(th_full * 180 / math.pi, q_full, color=ORANGE, lw=2.0, label="fully-catalytic")
    ax.set_xlabel(r"surface angle from stagnation point (deg)")
    ax.set_ylabel(r"$q_w$ (MW/m$^2$)")
    ax.set_title("wall heat flux bracket, Mach 20 N2 cylinder, isothermal 1000 K wall")
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

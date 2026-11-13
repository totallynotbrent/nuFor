# park5 cylinder gate figures on the shared house template: the
# stagnation-line profiles. regenerates from the run's vtk.
#
# the shock-layer fields themselves live in docs/img, plotted one panel
# per png by cases/fieldplots.py.
import sys
import os
import re
import math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figtemplate import figure, save, ORANGE, PURPLE, GRAY

import numpy as np


def read_fields(path):
    txt = open(path).read()
    blocks = re.findall(
        r'SCALARS (\S+) double 1\nLOOKUP_TABLE default\n((?:[-\d.eE+ ]+\n)+)',
        txt,
    )
    return {name: np.array([float(x) for x in body.split()]) for name, body in blocks}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    vtk = os.path.join(here, "cylinder-m20-park5.vtk")
    f = read_fields(vtk)
    nx, ny = 400, 200
    x0, x1, y0, y1 = -1.8, 5.0, 0.0, 3.4
    dx, dy = (x1 - x0) / nx, (y1 - y0) / ny

    # stagnation line (first row above the symmetry axis)
    j = 0
    xs, ts, xs_frac = [], [], []
    m_n2, m_n = 28.016e-3, 14.008e-3
    for i in range(nx):
        x = x0 + (i + 0.5) * dx
        if -1.5 < x < -0.9:
            k = i + j * nx
            y = f["y_N"][k]
            xs.append(x)
            ts.append(f["temperature"][k])
            xs_frac.append((y / m_n) / ((y / m_n) + (1 - y) / m_n2))
    fig2, (axt, axx) = figure(nrows=2, figsize=(7, 6))
    axt.plot(xs, ts, color=ORANGE)
    axt.axhspan(15000, 17300, color=GRAY, alpha=0.3)
    axt.set_ylabel("T (K)")
    axt.set_title("stagnation-line temperature (gray band: published one-T expectation)")
    axx.plot(xs, xs_frac, color=PURPLE)
    axx.set_xlabel("x (m)")
    axx.set_ylabel(r"$X_\mathrm{N}$ (mole fraction)")
    axx.set_title("atomic nitrogen along the stagnation line")
    save(fig2, "cyl-park5-stagline", out_dir=here)


if __name__ == "__main__":
    main()

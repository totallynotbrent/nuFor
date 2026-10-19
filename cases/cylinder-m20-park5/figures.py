# park5 cylinder gate figures on the shared house template: the
# shock-layer fields (temperature + atomic nitrogen) and the
# stagnation-line profiles. regenerates from the run's vtk.
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

    def grid(name):
        v = f[name]
        return np.array([v[i + j * nx] for j in range(ny) for i in range(nx)]).reshape(ny, nx)

    tg, yg = grid("temperature"), grid("y_N")
    body = tg <= 100.0
    tg_m = np.ma.masked_where(body, tg)
    yg_m = np.ma.masked_where(body, yg)
    xg = x0 + (np.arange(nx) + 0.5) * dx
    ycg = y0 + (np.arange(ny) + 0.5) * dy
    th = np.linspace(0, math.pi, 240)

    fig, axes = figure(ncols=2, figsize=(11, 2.8))
    im0 = axes[0].pcolormesh(xg, ycg, tg_m, cmap="inferno", vmax=17000)
    axes[0].plot(np.cos(th), np.sin(th), color="white", lw=1.2)
    axes[0].set_xlim(-1.8, 2.2)
    axes[0].set_ylim(0, 1.8)
    axes[0].set_aspect("equal")
    axes[0].set_title("temperature (K)")
    fig.colorbar(im0, ax=axes[0], shrink=0.85)
    im1 = axes[1].pcolormesh(xg, ycg, yg_m, cmap="viridis", vmax=0.5)
    axes[1].plot(np.cos(th), np.sin(th), color="white", lw=1.2)
    axes[1].set_xlim(-1.8, 2.2)
    axes[1].set_ylim(0, 1.8)
    axes[1].set_aspect("equal")
    axes[1].set_title(r"$Y_\mathrm{N}$ (mass fraction)")
    fig.colorbar(im1, ax=axes[1], shrink=0.85)
    fig.suptitle("Mach 20 nitrogen cylinder, park5 chemistry")
    fig.tight_layout()
    save(fig, "cyl-park5-fields", out_dir=here)

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

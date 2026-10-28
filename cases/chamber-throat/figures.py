"""the chamber-throat verification figures.

two figures on the house template:

- the settling: the throat's mach, pressure and temperature against the
  clock, and the chamber's feed against the throat's throughput, with
  cea's equilibrium throat state marked.
- the profile: the duct's axial pressure, temperature and mach at the
  settled state, with the wall's radius beside them.

the inputs are the run's interval snapshots and the throat history the
watcher writes beside them:

    throat-history.csv   step,t,throat_col,throat_area,throat_mdot,
                         throat_p,throat_T,throat_M,chamber_mdot,
                         chamber_p,chamber_T,chamber_M
    vtk-XXXXXX.vtk       one snapshot per interval

usage: figures.py [run-dir] [out-dir]
"""
import ast
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY  # noqa: E402
from chamber_fields import (  # noqa: E402
    read_case,
    read_history,
    read_vtk,
    columns,
    throat_column,
    wall_radius,
)

HERE = os.path.dirname(os.path.abspath(__file__))

# cea's equilibrium throat state for this chamber, and the closure's own
# frozen-composition isentrope at the chamber's pressure.
CEA_P = 4.05352e6
CEA_T = 3501.98
CLOSURE_RATIO = 0.56534


def settling_figure(history, out_dir):
    """the throat's state against the clock."""
    steps = [r["step"] for r in history]
    times = [r["t"] for r in history]
    fig, axes = figure(nrows=2, ncols=2, figsize=(9.5, 6.4))
    ax = axes[0][0]
    ax.plot(times, [r["throat_M"] for r in history], "o-", color=ORANGE, ms=3,
            label="throat mach")
    ax.axhline(1.0, color=GRAY, lw=1.0, label="mach 1")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("mach")
    ax.set_title("the throat chokes")
    ax.legend(loc="best", fontsize=9)

    ax = axes[0][1]
    ax.plot(times, [r["throat_p"] / 1e6 for r in history], "o-", color=PURPLE, ms=3,
            label="throat")
    ax.axhline(CEA_P / 1e6, color=RED, lw=1.0, label="cea throat")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("pressure (MPa)")
    ax.set_title("the throat's pressure")
    ax.legend(loc="best", fontsize=9)

    ax = axes[1][0]
    ax.plot(times, [r["throat_T"] for r in history], "o-", color=GREEN, ms=3,
            label="throat")
    ax.axhline(CEA_T, color=RED, lw=1.0, label="cea throat")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("temperature (K)")
    ax.set_title("the throat's temperature")
    ax.legend(loc="best", fontsize=9)

    ax = axes[1][1]
    ax.plot(times, [r["chamber_mdot"] for r in history], "o-", color=ORANGE, ms=3,
            label="chamber feed")
    ax.plot(times, [r["throat_mdot"] for r in history], "s-", color=PURPLE, ms=3,
            label="throat throughput")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("mass flow (kg/s)")
    ax.set_title("the feed against the throat's throughput")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    save(fig, "chamber-throat-settling", out_dir=out_dir)


def profile_figure(case, cols, out_dir):
    """the duct's axial state at the settled snapshot."""
    nx = case["nx"]
    ti, min_r = throat_column(case, cols)
    dx = (case["x1"] - case["x0"]) / nx
    xs = [case["x0"] + (i + 0.5) * dx for i in sorted(cols)]
    fig, axes = figure(nrows=5, figsize=(7.6, 10.2))
    ax = axes[0]
    radii = [wall_radius(case["verts"], x) for x in xs]
    ax.plot(xs, radii, color=GRAY, lw=1.4)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("radius (m)")
    ax.set_title("the duct's wall")

    for ax, key, label, title, color in (
        (axes[1], "p", "pressure (MPa)", "the axial pressure", PURPLE),
        (axes[2], "T", "temperature (K)", "the axial temperature", GREEN),
        (axes[3], "mach", "mach", "the axial mach", ORANGE),
    ):
        ys = []
        for i in sorted(cols):
            v = cols[i][key]
            ys.append(v / 1e6 if key == "p" else v)
        ax.plot(xs, ys, color=color, lw=1.4)
        ax.axvline(case["x0"] + (ti + 0.5) * dx, color=RED, lw=1.0, ls="--")
        ax.set_xlabel("x (m)")
        ax.set_ylabel(label)
        ax.set_title(title)
    axes[3].axhline(1.0, color=GRAY, lw=1.0, ls=":", label="mach 1")
    axes[1].axhline(CEA_P / 1e6, color=RED, lw=1.0, ls=":", label="cea throat")
    axes[2].axhline(CEA_T, color=RED, lw=1.0, ls=":", label="cea throat")
    axes[1].legend(loc="best", fontsize=9)
    axes[2].legend(loc="best", fontsize=9)
    axes[3].legend(loc="best", fontsize=9)

    ax = axes[4]
    ys = [cols[i]["m_dot"] for i in sorted(cols)]
    ax.plot(xs, ys, color=PURPLE, lw=1.4)
    ax.axvline(case["x0"] + (ti + 0.5) * dx, color=RED, lw=1.0, ls="--")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("mass flow (kg/s)")
    ax.set_title("the axial mass flow")
    axes[0].axvline(case["x0"] + (ti + 0.5) * dx, color=RED, lw=1.0, ls="--")
    axes[0].text(case["x0"] + (ti + 0.5) * dx, max(radii) * 0.92, " throat",
                 color=RED, fontsize=9)
    fig.tight_layout()
    save(fig, "chamber-throat-profile", out_dir=out_dir)


def main():
    run_dir = sys.argv[1] if len(sys.argv) > 1 else "run"
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        HERE, "..", "..", "docs", "numerics", "2d"
    )
    case = read_case(os.path.join(HERE, "case.toml"))
    history = read_history(os.path.join(run_dir, "throat-history.csv"))
    settling_figure(history, out_dir)
    snapshots = sorted(f for f in os.listdir(run_dir) if f.startswith("vtk-"))
    fields, nx, ny, origin, spacing = read_vtk(os.path.join(run_dir, snapshots[-1]))
    cols = columns(case, fields, nx, ny, origin, spacing)
    profile_figure(case, cols, out_dir)


if __name__ == "__main__":
    main()

"""the chamber-throat verification figures.

three figures on the house template:

- the settling: the sonic pressure, temperature and mass flow against
  the clock, with the marched gas's own isentrope and cea's equilibrium
  throat state marked.
- the profile: the duct's axial pressure, temperature and mach at the
  settled state, with the wall's radius beside them.
- the loss: the marched pressure against the isentrope the gas's own
  thermochemistry demands, which is where the marcher's error shows.

the inputs are the distilled record the case ships:

    settled/history.csv    step,sonic_col,sonic_x_m,sonic_p_pa,sonic_T_k,
                           sonic_mdot_kg_s,chamber_p_pa,chamber_T_k,
                           chamber_mdot_kg_s
    settled/profile.csv    x_m,area_m2,p_pa,T_k,mixed_mach,flux_mach,
                           velocity_m_s,m_dot_kg_s
    settled/isentrope.csv  x_m,isentropic_p_pa,marched_p_pa,deviation_percent

usage: figures.py [settled-dir] [out-dir]
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# cea's equilibrium throat state for this chamber
CEA_P = 4.05352e6
CEA_T = 3501.98


def read_rows(path):
    with open(path, newline="") as handle:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(row for row in handle if not row.startswith("#"))
        ]


def settling_figure(history, isentrope, out_dir):
    """the sonic state against the clock."""
    steps = [row["step"] for row in history]
    fig, axes = figure(nrows=2, ncols=2, figsize=(9.5, 6.4))

    ax = axes[0][0]
    ax.plot(steps, [row["sonic_p_pa"] / 1e6 for row in history], "o-", color=PURPLE, ms=3,
            label="the sonic point")
    ax.axhline(CEA_P / 1e6, color=RED, lw=1.0, ls=":", label="cea throat")
    if isentrope:
        ax.axhline(isentrope[0]["isentropic_p_pa"] / 1e6, color=GRAY, lw=1.2, ls="--",
                   label="the gas's own isentrope")
    ax.set_xlabel("step")
    ax.set_ylabel("pressure (MPa)")
    ax.set_title("the sonic pressure")
    ax.legend(loc="best", fontsize=9)

    ax = axes[0][1]
    ax.plot(steps, [row["sonic_T_k"] for row in history], "o-", color=GREEN, ms=3,
            label="the sonic point")
    ax.axhline(CEA_T, color=RED, lw=1.0, ls=":", label="cea throat")
    ax.set_xlabel("step")
    ax.set_ylabel("temperature (K)")
    ax.set_title("the sonic temperature")
    ax.legend(loc="best", fontsize=9)

    ax = axes[1][0]
    ax.plot(steps, [row["sonic_mdot_kg_s"] for row in history], "s-", color=ORANGE, ms=3,
            label="through the throat")
    ax.axhline(3.80, color=RED, lw=1.0, ls=":", label="cea")
    ax.set_xlabel("step")
    ax.set_ylabel("mass flow (kg/s)")
    ax.set_title("the throat's throughput")
    ax.legend(loc="best", fontsize=9)

    ax = axes[1][1]
    ax.plot(steps, [row["chamber_p_pa"] / 1e6 for row in history], "o-", color=PURPLE, ms=3,
            label="the chamber")
    ax.set_xlabel("step")
    ax.set_ylabel("pressure (MPa)")
    ax.set_title("the chamber's pressure")
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    save(fig, "chamber-throat-settling", out_dir=out_dir)


def profile_figure(profile, sonic_x, out_dir):
    """the duct's axial state at the settled snapshot."""
    xs = [row["x_m"] for row in profile]
    fig, axes = figure(nrows=5, figsize=(7.6, 10.2))

    ax = axes[0]
    areas = [row["area_m2"] for row in profile]
    ax.plot(xs, [math.sqrt(area / math.pi) for area in areas], color=GRAY, lw=1.4)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("equivalent radius (m)")
    ax.set_title("the duct's flow area")

    for ax, key, scale, label, title, color in (
        (axes[1], "p_pa", 1e-6, "pressure (MPa)", "the axial pressure", PURPLE),
        (axes[2], "T_k", 1.0, "temperature (K)", "the axial temperature", GREEN),
        (axes[3], "mixed_mach", 1.0, "mach", "the axial mach", ORANGE),
    ):
        ax.plot(xs, [row[key] * scale for row in profile], color=color, lw=1.4)
        ax.axvline(sonic_x, color=RED, lw=1.0, ls="--")
        ax.set_xlabel("x (m)")
        ax.set_ylabel(label)
        ax.set_title(title)
    axes[3].plot(xs, [row["flux_mach"] for row in profile], color=PURPLE, lw=1.0, ls=":",
                 label="the flux's mach")
    axes[3].axhline(1.0, color=GRAY, lw=1.0, ls=":", label="mach 1")
    axes[1].axhline(CEA_P / 1e6, color=RED, lw=1.0, ls=":", label="cea throat")
    axes[2].axhline(CEA_T, color=RED, lw=1.0, ls=":", label="cea throat")
    axes[1].legend(loc="best", fontsize=9)
    axes[2].legend(loc="best", fontsize=9)
    axes[3].legend(loc="best", fontsize=9)

    ax = axes[4]
    ax.plot(xs, [row["m_dot_kg_s"] for row in profile], color=PURPLE, lw=1.4)
    ax.axvline(sonic_x, color=RED, lw=1.0, ls="--")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("mass flow (kg/s)")
    ax.set_title("the axial mass flow")
    axes[0].axvline(sonic_x, color=RED, lw=1.0, ls="--")
    axes[0].text(sonic_x, max(math.sqrt(a / math.pi) for a in areas) * 0.92, " sonic",
                 color=RED, fontsize=9)
    fig.tight_layout()
    save(fig, "chamber-throat-profile", out_dir=out_dir)


def loss_figure(isentrope, out_dir):
    """the marched pressure against the isentrope the gas demands."""
    xs = [row["x_m"] * 1000.0 for row in isentrope]
    fig, axes = figure(nrows=2, figsize=(7.6, 5.6))

    ax = axes[0]
    ax.plot(xs, [row["isentropic_p_pa"] / 1e6 for row in isentrope], "-", color=GRAY, lw=1.4,
            label="the gas's own isentrope")
    ax.plot(xs, [row["marched_p_pa"] / 1e6 for row in isentrope], "-", color=PURPLE, lw=1.4,
            label="the march")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("pressure (MPa)")
    ax.set_title("the marched pressure against the isentrope")
    ax.legend(loc="best", fontsize=9)

    ax = axes[1]
    deviations = [row["deviation_percent"] for row in isentrope]
    ax.plot(xs, deviations, "-", color=ORANGE, lw=1.4)
    ax.axhline(0.0, color=GRAY, lw=1.0)
    ax.set_ylim(min(deviations) - 1.0, 1.0)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("percent off the isentrope")
    ax.set_title("how far the march sits off the isentrope")
    fig.tight_layout()
    save(fig, "chamber-throat-isentrope-loss", out_dir=out_dir)


def main():
    settled = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "settled")
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        HERE, "..", "..", "docs", "numerics", "2d"
    )
    history = read_rows(os.path.join(settled, "history.csv"))
    profile = read_rows(os.path.join(settled, "profile.csv"))
    isentrope = read_rows(os.path.join(settled, "isentrope.csv"))
    sonic_x = history[-1]["sonic_x_m"]
    settling_figure(history, isentrope, out_dir)
    profile_figure(profile, sonic_x, out_dir)
    loss_figure(isentrope, out_dir)


if __name__ == "__main__":
    main()

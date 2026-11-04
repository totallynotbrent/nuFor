#!/usr/bin/env python3
"""the ignition figures for the chemistry page.

these figures back the zero-dimensional ignition claim with the
numbers the case actually produced:

- the temperature history of the burnt parcel against cantera's own
  constant-volume reactor for the same parcel.
- the flash window on its own axes, where the two integrators are
  visibly separated and the separation is the discretisation error.
- the major species histories, cantera's drawn underneath as a grey
  underlay so the agreement is judged by the coloured curve.
- the temperature error through the whole burn, on a logarithmic
  axis, for the marched step and a coarser one.

usage:

    env -u PYTHONPATH /usr/bin/python3 ignition_figures.py \\
        <case-dir> <out-dir>

the case dir holds ignition-trace.csv, the marching output, and
cantera-reference.csv, the oracle.
"""

import bisect
import csv
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY  # noqa: E402

SPECIES = [
    ("CH4", PURPLE),
    ("O2", GREEN),
    ("H2O", "#2E86C1"),
    ("CO2", ORANGE),
    ("CO", RED),
]


def read_csv(path):
    """read a trace into column dictionaries."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        raise SystemExit(f"{path} has no rows")
    columns = {}
    for key in rows[0]:
        columns[key] = [float(row[key]) for row in rows]
    return columns


def interpolate(xs, ys, at):
    """the linear interpolation of ys at at, clamped at the ends."""
    if at <= xs[0]:
        return ys[0]
    if at >= xs[-1]:
        return ys[-1]
    i = bisect.bisect_left(xs, at)
    lo, hi = xs[i - 1], xs[i]
    span = hi - lo
    if span <= 0.0:
        return ys[i]
    weight = (at - lo) / span
    return ys[i - 1] * (1.0 - weight) + ys[i] * weight


def style(ax):
    """the house axes style."""
    ax.grid(True, alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)


def plot_temperature(history, ours, theirs, out_dir):
    """the whole burn, ours against cantera's."""
    fig, ax = figure(figsize=(7.2, 4.4))
    ax.plot(theirs["time_s"], theirs["temperature_K"],
            color=GRAY, linewidth=3.0, label="Cantera 3.2", zorder=1)
    ax.plot(ours["time_s"], ours["temperature_K"],
            color=ORANGE, linewidth=1.5, label="nuFor, marched", zorder=3)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("temperature (K)")
    ax.set_title("Stoichiometric CH$_4$/air at 2000 K, rigid vessel")
    ax.legend(loc="center right")
    style(ax)
    fig.tight_layout()
    return save(fig, "gri30-ignition-temperature", out_dir)


def plot_flash_window(history, ours, theirs, out_dir):
    """the flash window, where the discretisation error is visible."""
    fig, ax = figure(figsize=(7.2, 4.4))
    ax.plot(theirs["time_s"], theirs["temperature_K"],
            color=GRAY, linewidth=3.2, label="Cantera 3.2", zorder=1)
    ax.plot(ours["time_s"], ours["temperature_K"],
            color=ORANGE, linewidth=1.5, label="nuFor, marched", zorder=3)
    ax.set_xlim(0.0, 1.2e-4)
    ax.set_ylim(1900.0, 3060.0)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("temperature (K)")
    ax.set_title("Flash window, enlarged")
    ax.legend(loc="lower right")
    style(ax)
    fig.tight_layout()
    return save(fig, "gri30-ignition-flash-window", out_dir)


def plot_species(history, ours, theirs, out_dir):
    """the major species histories, with cantera's as a grey underlay."""
    fig, ax = figure(figsize=(7.2, 4.4))
    for name, color in SPECIES:
        key = f"massfrac_{name}"
        if key not in ours or key not in theirs:
            continue
        ax.plot(theirs["time_s"], theirs[key],
                color="#cfcfcf", linewidth=3.0, zorder=1)
        ax.plot(ours["time_s"], ours[key],
                color=color, linewidth=1.4, label=name, zorder=3)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("mass fraction")
    ax.set_title("Major species, grey underlay is Cantera")
    ax.set_ylim(-0.006, 0.16)
    ax.legend(loc="upper center", ncol=5)
    style(ax)
    fig.tight_layout()
    return save(fig, "gri30-ignition-species", out_dir)


def plot_error(history, ours, theirs, coarse, out_dir):
    """the temperature error through the burn, for two step sizes."""
    fig, ax = figure(figsize=(7.2, 4.4))
    samples = theirs["time_s"]
    curves = []

    ours_times = ours["time_s"]
    ours_temps = ours["temperature_K"]
    error = [
        abs(interpolate(ours_times, ours_temps, t) - theirs["temperature_K"][i])
        for i, t in enumerate(samples)
    ]
    median = sorted(error)[len(error) // 2]
    ax.semilogy(samples, [max(e, 1.0e-3) for e in error],
                color=ORANGE, linewidth=1.4,
                label=f"step 0.5 us, median {median:.4f} K")
    curves.append(("step 0.5 us", median, max(error)))

    if coarse is not None:
        coarse_times = coarse["time_s"]
        coarse_temps = coarse["temperature_K"]
        coarse_error = [
            abs(interpolate(coarse_times, coarse_temps, t) - theirs["temperature_K"][i])
            for i, t in enumerate(samples)
        ]
        coarse_median = sorted(coarse_error)[len(coarse_error) // 2]
        ax.semilogy(samples, [max(e, 1.0e-3) for e in coarse_error],
                    color=PURPLE, linewidth=1.1, linestyle="--",
                    label=f"step 2 us, median {coarse_median:.4f} K")
        curves.append(("step 2 us", coarse_median, max(coarse_error)))

    # the flash dominates the error; split the axis at the point where
    # the parcel has finished burning so both halves can be read.
    settled = [(t, v) for t, v in zip(samples, error) if t > 2.0e-4]
    if settled:
        post = sorted(v for _, v in settled)
        post_median = post[len(post) // 2]
        ax.axvspan(0.0, 2.0e-4, color="#eeeeee", zorder=0)
        ax.annotate(
            f"flash window: transient error\nsettled after 0.2 ms, median {post_median:.3f} K",
            xy=(2.05e-4, 3.0),
            xytext=(5.0e-4, 12.0),
            fontsize=9,
            color="#333333",
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#999999", linewidth=0.8),
            arrowprops=dict(arrowstyle="-", color="#888888", linewidth=0.8),
        )
    ax.set_xlabel("time (s)")
    ax.set_ylabel("|nuFor $-$ Cantera| temperature (K)")
    ax.set_title("Temperature error through the flash")
    ax.legend(loc="upper right")
    style(ax)
    fig.tight_layout()

    history.extend(curves)
    return save(fig, "gri30-ignition-error", out_dir)


def main(argv):
    """draw every figure and report the numbers the page quotes."""
    if len(argv) != 3:
        print(__doc__)
        return 2
    case_dir = argv[1]
    out_dir = argv[2]

    ours = read_csv(os.path.join(case_dir, "ignition-trace.csv"))
    theirs = read_csv(os.path.join(case_dir, "cantera-reference.csv"))
    coarse_path = os.path.join(case_dir, "ignition-trace-coarse.csv")
    coarse = read_csv(coarse_path) if os.path.exists(coarse_path) else None

    history = []
    plot_temperature(history, ours, theirs, out_dir)
    plot_flash_window(history, ours, theirs, out_dir)
    plot_species(history, ours, theirs, out_dir)
    plot_error(history, ours, theirs, coarse, out_dir)

    end = min(ours["time_s"][-1], theirs["time_s"][-1])
    our_temp = interpolate(ours["time_s"], ours["temperature_K"], end)
    ct_temp = interpolate(theirs["time_s"], theirs["temperature_K"], end)
    our_press = interpolate(ours["time_s"], ours["pressure_Pa"], end)
    ct_press = interpolate(theirs["time_s"], theirs["pressure_Pa"], end)
    print(f"final temperature: nuFor {our_temp:.4f} K, Cantera {ct_temp:.4f} K,"
          f" difference {our_temp - ct_temp:+.4f} K")
    print(f"final pressure: nuFor {our_press:.3f} Pa, Cantera {ct_press:.3f} Pa,"
          f" difference {our_press - ct_press:+.3f} Pa")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

"""the engine-wall bracket's figures.

the bracket isolates the wall's contribution to the chamber's thermal
balance by changing ONLY the wall's treatment and marching every limb
to the same clock:

- the steadiness: the residual, the density maximum and the peak
  temperature against the clock for every limb, so the divergence that
  used to kill a cooled wall is visible as a curve.
- the duty: the power each wall treatment draws out of the gas, split
  into its conduction, radiation and film shares.
- the stations: the duty distribution along the chamber's axis for
  every limb.

the inputs are the case's own recorded histories:

    <leg>/run-log.csv      step,t,...

usage: figures.py [runs-dir] [out-dir]

the figure library lives in the system's python, not in an agent's
virtual environment, so run this with the system interpreter and a clean
PYTHONPATH:

    env -u PYTHONPATH /usr/bin/python3 figures.py <runs-dir> <out-dir>
"""
import csv
import json
import os
import sys

import numpy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# the limbs of the bracket, in the order they are compared. the first
# field is the short name the figures use, the second the name on the
# plot, and the third the colour.
#
# the runs live in directories named either after the short name or
# after the case's tag, so the resolver tries both.
LIMBS = (
    ("adiabatic", "adiabatic", GRAY),
    ("radiating", "radiating", RED),
    ("filmcooled", "film cooled", ORANGE),
    ("regen", "regenerative", PURPLE),
    ("conjugate", "conjugate", GREEN),
)

# the case tag each short name travels under, when it travels under one.
LEG_TAGS = {
    "adiabatic": "wall-adiabatic",
    "radiating": "wall-radiating",
    "filmcooled": "wall-filmcooled",
    "regen": "wall-regenerative",
    "conjugate": "wall-conjugate",
}


def resolve_run_dir(runs_dir, name):
    """the run's directory, whichever of its names it answers to."""
    for candidate in (name, LEG_TAGS.get(name)):
        if not candidate:
            continue
        path = os.path.join(runs_dir, candidate)
        if os.path.isdir(path):
            return path
    return os.path.join(runs_dir, name)



def panels(axes):
    """the figure's axes as a flat list, whatever shape matplotlib gave."""
    return list(numpy.atleast_1d(axes).reshape(-1))


def one_panel(axes):
    """a single axes, whether matplotlib gave it bare or wrapped."""
    if hasattr(axes, "set_xlabel"):
        return axes
    return panels(axes)[0]

def read_samples(run_dir):
    """the run's recorded samples, oldest first."""
    path = os.path.join(run_dir, "progress.json")
    if not os.path.exists(path):
        return []
    with open(path) as handle:
        data = json.load(handle)
    samples = data.get("samples") or []
    if not samples and isinstance(data.get("sample"), dict):
        samples = [data["sample"]]
    return samples


def read_log(run_dir):
    """the run's summary lines, joined."""
    path = os.path.join(run_dir, "run-log.txt")
    if not os.path.exists(path):
        return ""
    with open(path) as handle:
        return handle.read()


def parse_duty(log):
    """the wall duty's total and its three shares, in MW."""
    duty = {"total": 0.0, "conduction": 0.0, "radiation": 0.0, "film": 0.0}
    for line in log.splitlines():
        stripped = line.strip()
        if stripped.startswith("wall duty:"):
            duty["total"] = float(stripped.split(":")[1].split()[0])
        elif stripped.startswith("(conduction"):
            parts = stripped.replace("(", "").replace(")", "").split(",")
            for part in parts:
                key, _, value = part.strip().partition(" ")
                duty[key] = float(value.split()[0])
    return duty


def parse_surface_temperature(log):
    """the temperature the gas met at the wall, in K, if it was reported."""
    for line in log.splitlines():
        stripped = line.strip()
        if stripped.startswith("wall surface:"):
            parts = stripped.replace(",", " ").split()
            try:
                return float(parts[2])
            except (IndexError, ValueError):
                return None
    return None


def parse_stations(log):
    """the duty's distribution along the axis, as (x0, x1, mw) triples."""
    stations = []
    in_block = False
    for line in log.splitlines():
        stripped = line.strip()
        if stripped.startswith("wall duty by station"):
            in_block = True
            continue
        if not in_block:
            continue
        if not stripped or stripped.startswith(("wall ", "(", "aero")):
            break
        parts = stripped.split()
        # the station lines read: x <lo>.. <hi> <mw> perimeter <metres> m
        if len(parts) < 4 or parts[0] != "x":
            break
        lo = parts[1].rstrip("..")
        hi = parts[2]
        try:
            stations.append((float(lo), float(hi), float(parts[3]), float(parts[5])))
        except ValueError:
            break
    return stations


def steadiness_figure(runs, out_dir):
    """the residual, the density maximum and the peak temperature."""
    fig, axes = figure(nrows=1, ncols=3, figsize=(11.6, 3.9))
    axes = panels(axes)

    for name, label, color in LIMBS:
        samples = runs[name]["samples"]
        if not samples:
            continue
        steps = [row["step"] for row in samples]
        times = [row["t"] for row in samples]
        axes[0].plot(times, [row["residual"] for row in samples], "-", color=color,
                     lw=1.3, label=label)
        axes[1].plot(times, [row["rho_max"] for row in samples], "-", color=color,
                     lw=1.3, label=label)
        axes[2].plot(times, [row["t_max"] for row in samples], "-", color=color,
                     lw=1.3, label=label)

    axes[0].axhline(5.0, color=RED, lw=1.0, ls=":", label="the divergence mark")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("time (s)")
    axes[0].set_ylabel("residual")
    axes[0].set_title("the residual against the clock")
    axes[0].legend(loc="best", fontsize=8)

    axes[1].set_xlabel("time (s)")
    axes[1].set_ylabel("density maximum (kg/m3)")
    axes[1].set_title("the density maximum")
    axes[1].legend(loc="best", fontsize=8)

    axes[2].set_xlabel("time (s)")
    axes[2].set_ylabel("temperature (K)")
    axes[2].set_title("the peak temperature")
    axes[2].legend(loc="best", fontsize=8)
    fig.tight_layout()
    save(fig, "engine-wall-bracket-steadiness", out_dir=out_dir)


def duty_figure(runs, out_dir):
    """the power each wall treatment draws out of the gas."""
    fig, axes = figure(nrows=1, ncols=2, figsize=(9.5, 3.9))
    axes = panels(axes)

    drawn = [(name, label, color) for name, label, color in LIMBS
             if runs[name]["duty"]["total"] > 0.0]
    names = [label for _, label, _ in drawn]
    colors = [color for _, _, color in drawn]
    widths = 0.34
    positions = list(range(len(names)))

    bottom = [0.0] * len(names)
    for key, label, color in (
        ("conduction", "conduction", GRAY),
        ("radiation", "radiation", ORANGE),
        ("film", "film", PURPLE),
    ):
        values = [runs[name]["duty"][key] for name, _, _ in drawn]
        axes[0].bar([p + widths / 2.0 for p in positions], values, width=widths,
                    bottom=bottom, color=color, label=label)
        bottom = [b + v for b, v in zip(bottom, values)]

    axes[0].set_xticks([p + widths / 2.0 for p in positions])
    axes[0].set_xticklabels(names, fontsize=8)
    axes[0].set_ylabel("power (MW)")
    axes[0].set_title("the duty each wall draws from the gas")
    axes[0].legend(loc="best", fontsize=9)

    axes[1].bar(names, [runs[name]["duty"]["total"] for name, _, _ in drawn],
                color=colors, alpha=0.85)
    axes[1].set_ylabel("power (MW)")
    axes[1].set_title("the total duty")
    axes[1].tick_params(axis="x", labelrotation=30)
    for tick in axes[1].get_xticklabels():
        tick.set_fontsize(8)
    fig.tight_layout()
    save(fig, "engine-wall-bracket-duty", out_dir=out_dir)


def surface_figure(runs, out_dir):
    """the temperature the gas actually met at each wall."""
    fig, axes = figure(nrows=1, ncols=1, figsize=(7.6, 4.2))
    axes = one_panel(axes)

    drawn = []
    for name, label, color in LIMBS:
        temperature = runs[name]["surface_k"]
        if temperature is None:
            continue
        drawn.append((name, label, color, temperature))

    if drawn:
        axes.bar(
            [label for _, label, _, _ in drawn],
            [temperature for _, _, _, temperature in drawn],
            color=[color for _, _, color, _ in drawn],
            alpha=0.85,
        )
        axes.axhline(400.0, color=GRAY, lw=1.0, ls=":",
                     label="the plain wall's stated temperature")
        for position, (_, _, _, temperature) in enumerate(drawn):
            axes.annotate("{:.0f}".format(temperature),
                          xy=(position, temperature),
                          xytext=(0, 4),
                          textcoords="offset points",
                          ha="center",
                          fontsize=8)
    axes.set_ylabel("temperature (K)")
    axes.set_title("the surface temperature the gas met")
    axes.tick_params(axis="x", labelrotation=30)
    for tick in axes.get_xticklabels():
        tick.set_fontsize(8)
    axes.legend(loc="best", fontsize=9)
    fig.tight_layout()
    save(fig, "engine-wall-bracket-surface", out_dir=out_dir)


def stations_figure(runs, out_dir):
    """the duty's distribution along the chamber's axis."""
    fig, axes = figure(nrows=1, ncols=1, figsize=(7.6, 4.2))
    axes = one_panel(axes)

    for name, label, color in LIMBS:
        stations = runs[name]["stations"]
        if not stations:
            continue
        centres = [0.5 * (lo + hi) for lo, hi, _, _ in stations]
        powers = [mw for _, _, mw, _ in stations]
        axes.plot(centres, powers, "o-", color=color, ms=4, lw=1.3, label=label)

    axes.set_xlabel("axial station (m)")
    axes.set_ylabel("duty (MW per station)")
    axes.set_title("where each wall treatment takes the duty")
    axes.legend(loc="best", fontsize=9)
    fig.tight_layout()
    save(fig, "engine-wall-bracket-stations", out_dir=out_dir)


def main():
    runs_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "runs")
    out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        HERE, "..", "..", "docs", "numerics", "2d"
    )

    runs = {}
    for name, _, _ in LIMBS:
        run_dir = resolve_run_dir(runs_dir, name)
        log = read_log(run_dir)
        runs[name] = {
            "samples": read_samples(run_dir),
            "log": log,
            "duty": parse_duty(log),
            "stations": parse_stations(log),
            "surface_k": parse_surface_temperature(log),
        }

    steadiness_figure(runs, out_dir)
    duty_figure(runs, out_dir)
    stations_figure(runs, out_dir)
    surface_figure(runs, out_dir)


if __name__ == "__main__":
    main()

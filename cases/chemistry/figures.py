"""the GRI-Mech 3.0 verification figures.

these figures back the mechanics page with the crate's own numbers.
the data comes from the rust example, which reads the mechanism's
tables and the checked-in cantera fixture:

    cargo run --release -p nufor-core --example gri30_verification \
        > verification.csv

three panels make the point:

- the net rates of progress, ours against cantera's, for every
  reaction at every evaluated mixture. the diagonal is exact
  agreement and the points sit on it.
- the troe broadening factors against the reduced pressure, for a
  representative falloff reaction at three temperatures. the sharp
  knees are where the crossover bends.
- the third-body bath for one mixture, split by collider kind, showing
  that a named collider draws only itself while a generic one draws
  every species weighted by its efficiency.

usage:

    env -u PYTHONPATH /usr/bin/python3 figures.py <csv> <out-dir>
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY  # noqa: E402


def read_tables(path):
    """read the example's csv into named tables."""
    tables = {}
    section = None
    header = None
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if not row:
                continue
            if row[0].startswith("#"):
                section = row[0].lstrip("# ").strip()
                tables[section] = []
                header = None
                continue
            if section is None:
                # the first table precedes any section marker
                section = "net rates of progress"
                tables[section] = []
                header = row
                continue
            if header is None:
                header = row
                continue
            tables[section].append(dict(zip(header, row)))
    return tables


def plot_net_rates(table, out_dir):
    """our net rates against cantera's, one point per reaction.

    the rates span nine orders of magnitude and many of them vanish,
    so the panel is drawn in the magnitude plane on logarithmic axes.
    a point on the diagonal means the magnitudes agree exactly; the
    sign is checked separately by the agreement figure, which ranks
    every reaction's relative difference.
    """
    points = []
    for row in table:
        mine = abs(float(row["ours"]))
        cantera = abs(float(row["cantera"]))
        if mine <= 0.0 or cantera <= 0.0:
            continue
        points.append((mine, cantera))

    fig, raw_axes = figure(figsize=(7.6, 5.6))
    ax = raw_axes
    lower = min(min(point) for point in points)
    upper = max(max(point) for point in points)
    ax.plot([lower, upper], [lower, upper], color=GRAY, linewidth=1.0,
            linestyle="--", label="exact agreement")
    ax.scatter([point[1] for point in points], [point[0] for point in points],
               s=13, color=ORANGE, alpha=0.55, linewidths=0,
               label=f"{len(points)} reactions")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("cantera net rate magnitude")
    ax.set_ylabel("nuFor net rate magnitude")
    ax.set_title("GRI-Mech 3.0 net rates: ours against cantera")
    ax.legend(loc="upper left")
    return save(fig, "gri30-net-rates-vs-cantera", out_dir)


def plot_relative_error(table, out_dir):
    """the relative disagreement for every reaction, per mixture."""
    mixtures = {}
    for row in table:
        mine = float(row["ours"])
        cantera = float(row["cantera"])
        if mine == 0.0 and cantera == 0.0:
            continue
        relative = abs(mine - cantera) / max(abs(cantera), 1e-300)
        mixtures.setdefault(int(row["mixture"]), []).append(relative)

    fig, raw_axes = figure(figsize=(7, 4.2))
    ax = raw_axes
    names = []
    colours = [ORANGE, PURPLE, GREEN, RED]
    for position, mixture in enumerate(sorted(mixtures)):
        values = sorted(mixtures[mixture])
        names.append(str(mixture))
        ax.semilogy(range(len(values)), [max(v, 1e-18) for v in values],
                    color=colours[position % len(colours)], linewidth=1.2,
                    marker="o", markersize=2.5, label=f"mixture {mixture}")
    ax.set_xlabel("reaction, ranked within its mixture")
    ax.set_ylabel("relative difference")
    ax.set_title("Agreement with cantera across every reaction and mixture")
    ax.set_ylim(1e-18, 1e-9)
    ax.legend(loc="lower right", fontsize=9)
    return save(fig, "gri30-cantera-agreement", out_dir)


def plot_broadening(table, out_dir):
    """the troe broadening factor against reduced pressure.

    the mechanism mixes two kinds of falloff reaction. a plain one
    blends its two limits smoothly, so its broadening factor is one by
    definition. a troe one carries four parameters and its factor bends
    sharply away from the crossover. both are plotted, because the
    contrast is the point: the troe curve departs from unity where the
    plain one cannot.
    """
    groups = {}
    kinds = {}
    for row in table:
        key = (int(row["reaction"]), float(row["temperature_K"]))
        groups.setdefault(key, []).append(
            (float(row["reduced_pressure"]), float(row["broadening_factor"])))
        kinds[key[0]] = row["kind"]

    troe_ids = sorted({rid for rid in kinds if kinds[rid] == "troe"})
    plain_ids = sorted({rid for rid in kinds if kinds[rid] == "plain"})
    chosen_troe = troe_ids[0] if troe_ids else None
    chosen_plain = plain_ids[0] if plain_ids else None

    fig, raw_axes = figure(figsize=(7.6, 5.0))
    ax = raw_axes
    colours = [ORANGE, PURPLE, GREEN]
    temperatures = sorted({temp for _, temp in groups})

    for position, temperature in enumerate(temperatures):
        points = groups.get((chosen_troe, temperature))
        if not points:
            continue
        points.sort()
        ax.loglog([p[0] for p in points], [p[1] for p in points],
                  color=colours[position % len(colours)], linewidth=1.7,
                  marker="o", markersize=3.2,
                  label=f"troe reaction {chosen_troe}, {temperature:.0f} K")

    if chosen_plain is not None:
        for position, temperature in enumerate(temperatures):
            points = groups.get((chosen_plain, temperature))
            if not points:
                continue
            points.sort()
            ax.loglog([p[0] for p in points], [p[1] for p in points],
                      color=GRAY, linewidth=1.0, linestyle=":",
                      label=f"plain reaction {chosen_plain} (unity)")

    ax.axhline(1.0, color=RED, linewidth=0.9, linestyle="--", alpha=0.6)
    ax.set_xlabel("reduced pressure")
    ax.set_ylabel("broadening factor")
    ax.set_title("Falloff broadening: troe bending away from the crossover")
    ax.legend(fontsize=8, loc="best")
    return save(fig, "gri30-troe-broadening", out_dir)


def plot_bath(table, out_dir):
    """the third-body bath, one reaction pair at a time.

    the two cases are picked so the contrast is direct: a generic-
    collider reaction and a named-collider reaction taken from the same
    mixture. each panel shows the share of the bath that each species
    brings, so the generic case spreads across many species and the
    named case concentrates on one. the shares make the panels
    comparable however different the absolute totals are.
    """
    generic = {}
    named = {}
    for row in table:
        bucket = named if row["kind"] == "named" else generic
        bucket.setdefault(int(row["reaction"]), []).append(
            (float(row["density_kmol_m3"]),
             float(row["efficiency"]),
             float(row["contribution"]),
             int(row["species"])))

    def pick(bucket):
        richest = max(bucket.items(),
                      key=lambda item: sum(entry[2] for entry in item[1]))
        return richest[0], sorted(richest[1], key=lambda entry: -entry[2])

    generic_id, generic_entries = pick(generic)
    named_id, named_entries = pick(named)

    fig, raw_axes = figure(nrows=1, ncols=2, figsize=(12.0, 5.2))
    axes = list(raw_axes)
    shades = plt.cm.Oranges([0.3 + 0.65 * index / 20.0 for index in range(21)])
    for ax, entries, identifier, title in (
        (axes[0], generic_entries, generic_id,
         "generic collider: every species joins"),
        (axes[1], named_entries, named_id,
         "named collider: only the named species joins"),
    ):
        total = sum(entry[2] for entry in entries)
        shares = [entry[2] / total * 100.0 if total > 0 else 0.0
                  for entry in entries]
        contributors = sum(1 for share in shares if share > 0.0)
        ax.bar(range(len(shares)), shares, color=shades[:len(shares)],
               edgecolor="white", linewidth=0.5)
        ax.set_xticks(range(len(shares)))
        ax.set_xticklabels([str(entry[3]) for entry in entries],
                           fontsize=7, rotation=90)
        ax.set_xlabel("species index")
        ax.set_ylabel("share of the bath, percent")
        ax.set_title(f"reaction {identifier}: {title}")
        ax.text(0.02, 0.94,
                f"{contributors} species contribute",
                transform=ax.transAxes, va="top", fontsize=10)
        ax.set_ylim(0.0, 105.0)
    fig.tight_layout()
    return save(fig, "gri30-third-body-bath", out_dir)


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    path = argv[1]
    out_dir = argv[2]
    tables = read_tables(path)
    plot_net_rates(tables["net rates of progress"], out_dir)
    plot_relative_error(tables["net rates of progress"], out_dir)
    plot_broadening(tables["broadening factors of the falloff reactions"], out_dir)
    plot_bath(tables["third-body bath contribution per species"], out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

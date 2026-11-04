"""the figures that show what the bracket distinguishes, and what it cannot.

two claims are argued here, each from the gate's own measurement path,
so a picture and the gate's verdict cannot drift apart.

the first figure argues that the wall treatments DISCRIMINATE, and that
the duty the solver reports and the cooling the gas shows are the same
thing seen twice. the wall's claim, converted into the temperature drop
it implies for the chamber, is set beside the drop the gas actually
registers at the wall. the second panel shows where that effect lives:
almost all of it is in the gas against the wall, and very little of it
reaches the chamber's bulk. the third panel carries the measured
choked-mass-flow departure with its uncertainty.

the second figure argues the opposite conclusion honestly. the chamber
does ring, the legs DO agree on its frequency and phase, and yet the
fitted mode cannot reconstruct any leg's own waveform, so it cannot be
removed and the absolute levels cannot be trusted on their own. this is
the evidence that overturned the earlier periodic-orbit story.

usage: effect_figures.py <runs-dir> <out-dir>

the figure library lives in the system's python, not in an agent's
virtual environment, so run this with the system interpreter and a
clean PYTHONPATH:

    env -u PYTHONPATH /usr/bin/python3 effect_figures.py \
        <runs-dir> <out-dir>
"""
import argparse
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import numpy
from matplotlib import pyplot

from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY  # noqa: E402

import cross_leg_settlement as cls
import gates

LABELS = {
    "wall-adiabatic": "adiabatic",
    "wall-radiating": "radiating",
    "wall-filmcooled": "film cool",
    "wall-regenerative": "regenerative",
    "wall-conjugate": "conjugate",
}

COLORS = {
    "wall-adiabatic": GRAY,
    "wall-radiating": RED,
    "wall-filmcooled": ORANGE,
    "wall-regenerative": PURPLE,
    "wall-conjugate": GREEN,
}

# the line styles let two curves be told apart when they overlap.
STYLES = {
    "wall-adiabatic": "-",
    "wall-radiating": "--",
    "wall-filmcooled": "-.",
    "wall-regenerative": ":",
    "wall-conjugate": (0, (5, 2, 1, 2)),
}


def frame(ax, title, xlabel, ylabel):
    """the house panel furniture."""
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)


def gather(runs_dir):
    """the gate's own measurement of every leg."""
    rows = []
    for tag, wall_kind in gates.BRACKET:
        info = gates.measure_one_leg(tag, runs_dir)
        info["specified_wall_k"] = wall_kind
        info["progress_path"] = os.path.join(runs_dir, tag, "progress.json")
        rows.append(info)
    return rows


def field_times(rows):
    """the marched time of each saved field, in milliseconds.

    the saved fields are numbered by the step at which they were
    written, so the time comes from the leg's own progress record,
    sampled densely and matched to the nearest saved step. the axis is
    then a real time rather than a field counter.
    """
    times = {}
    for row in rows:
        path = row.get("progress_path")
        samples = []
        if path and os.path.exists(path):
            with open(path) as handle:
                document = json.load(handle)
            samples = [
                (float(s["step"]), float(s["t"]))
                for s in document.get("samples") or []
                if s.get("step") is not None and s.get("t") is not None
            ]
        samples.sort()
        for observation in row.get("observations") or []:
            index = cls.field_index(observation)
            if index is None or not samples:
                continue
            stamp = min(samples, key=lambda item: abs(item[0] - index))[1]
            times.setdefault(index, stamp)
    return times


def paired_stats(entries, window_set):
    """the mean departure and its standard error over the judged window."""
    chosen = [e for e in entries if e["index"] in window_set]
    if not chosen:
        return None
    if len(chosen) < 2:
        return {"delta": chosen[0]["value"], "sem": 0.0,
                "level": chosen[0]["level"], "n": len(chosen)}
    values = [e["value"] for e in chosen]
    mean = statistics.fmean(values)
    variance = sum((v - mean) ** 2 for v in values) / (len(values) - 1)
    sem = variance ** 0.5 / len(values) ** 0.5
    level = statistics.fmean([e["level"] for e in chosen])
    return {"delta": mean, "sem": sem, "level": level, "n": len(chosen)}


def effect_figure(rows, result, out_dir):
    """what the wall did to the gas, and where that effect lives."""
    window_fields = result.get("window_fields") or []
    window_set = set(window_fields)
    times = field_times(rows)
    by_field = cls.series_by_field(
        rows, ("mass_flow", "chamber_p0", "chamber_t0", "near_wall_t"))
    cooled = [row["tag"] for row in rows if row["tag"] != cls.REFERENCE_TAG]
    dep = cls.departures_from(
        by_field, cls.REFERENCE_TAG, cooled,
        ("mass_flow", "chamber_p0", "chamber_t0", "near_wall_t"))

    fig, axes = figure(nrows=2, ncols=2, figsize=(9.4, 7.2))
    fig.suptitle(
        "the reported duty and the gas's own cooling, and where that "
        "cooling sits")

    # the duty the wall reports, as a chamber cooling, beside the cooling
    # the gas registers against the wall.
    ax = axes[0][0]
    claimed = []
    measured = []
    labels2 = []
    for tag in cooled:
        leg = next(row for row in rows if row["tag"] == tag)
        reference_mass_flow = statistics.fmean(
            list(by_field["mass_flow"][cls.REFERENCE_TAG].values()))
        reference_t0 = statistics.fmean(
            list(by_field["chamber_t0"][cls.REFERENCE_TAG].values()))
        scale = gates.predict_duty_temperature_drop(
            {"mass_flow": reference_mass_flow, "chamber_t0": reference_t0},
            leg.get("duty_mw"))
        near = paired_stats(dep[tag]["near_wall_t"], window_set)
        if scale is None or near is None:
            continue
        claimed.append(abs(scale["scale_dt0_k"]))
        measured.append(abs(near["delta"]))
        labels2.append(LABELS[tag])
    positions = numpy.arange(len(labels2))
    width = 0.36
    ax.bar(positions - width / 2.0, claimed, width=width, color=PURPLE,
           label="implied by the duty")
    ax.bar(positions + width / 2.0, measured, width=width, color=ORANGE,
           label="measured at the wall")
    frame(ax, "the implied duty against the measured cooling", "",
          "$\\Delta T$ against the adiabat (K)")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels2, fontsize=9)
    ax.legend(fontsize=9)
    for index, (left, right) in enumerate(zip(claimed, measured)):
        ax.text(index - width / 2.0, left + 1.0, "%.0f" % left, ha="center",
                fontsize=8)
        ax.text(index + width / 2.0, right + 1.0, "%.0f" % right, ha="center",
                fontsize=8)

    # where the effect lives: the near-wall gas against the chamber bulk.
    ax = axes[0][1]
    near_depths = []
    bulk_depths = []
    labels3 = []
    for tag in cooled:
        near = paired_stats(dep[tag]["near_wall_t"], window_set)
        bulk = paired_stats(dep[tag]["chamber_t0"], window_set)
        if near is None or bulk is None:
            continue
        near_depths.append(abs(near["delta"]))
        bulk_depths.append(abs(bulk["delta"]))
        labels3.append(LABELS[tag])
    positions = numpy.arange(len(labels3))
    width = 0.36
    ax.bar(positions - width / 2.0, near_depths, width=width, color=ORANGE,
           label="at the wall")
    ax.bar(positions + width / 2.0, bulk_depths, width=width, color=GREEN,
           label="in the chamber bulk")
    frame(ax, "the near-wall cooling against the chamber bulk", "",
          "$\\Delta T$ against the adiabat (K)")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels3, fontsize=9)
    ax.legend(fontsize=9)
    for index, (near_value, bulk_value) in enumerate(
            zip(near_depths, bulk_depths)):
        ax.text(index - width / 2.0, near_value + 1.0, "%.1f" % near_value,
                ha="center", fontsize=8)
        ax.text(index + width / 2.0, bulk_value + 1.0, "%.2f" % bulk_value,
                ha="center", fontsize=8)

    # the measured choked mass flow departure, with its uncertainty.
    ax = axes[1][0]
    departures = []
    errors = []
    labels4 = []
    for tag in cooled:
        stats = paired_stats(dep[tag]["mass_flow"], window_set)
        if stats is None:
            continue
        departures.append(100.0 * stats["delta"] / stats["level"])
        errors.append(100.0 * stats["sem"] / stats["level"])
        labels4.append(LABELS[tag])
    positions = numpy.arange(len(labels4))
    ax.bar(positions, departures, yerr=errors, capsize=4,
           color=[COLORS[t] for t in cooled][: len(departures)],
           error_kw={"ecolor": "black", "elinewidth": 1.0})
    frame(ax, "the chilled mass flow, against the adiabat", "",
          "change in mass flow (%)")
    ax.set_xticks(positions)
    ax.set_xticklabels(labels4, fontsize=9)
    ax.axhline(0.0, color="black", linewidth=0.8)
    for index, value in enumerate(departures):
        ax.text(index, value - 0.03, "%.3f" % value, ha="center", va="top",
                fontsize=8)

    # the near-wall departure through the march, on a real time axis.
    ax = axes[1][1]
    window_stamps = [times.get(index) for index in window_fields]
    stamps = [stamp for stamp in window_stamps if stamp is not None]
    if stamps:
        ax.axvspan(1000.0 * min(stamps), 1000.0 * max(stamps),
                   color="#E8E8E8", zorder=0)
    for tag in cooled:
        entries = dep[tag]["near_wall_t"]
        xs = []
        ys = []
        for entry in entries:
            stamp = times.get(entry["index"])
            if stamp is None:
                continue
            xs.append(1000.0 * stamp)
            ys.append(entry["value"])
        ax.plot(xs, ys, linestyle=STYLES[tag], color=COLORS[tag],
                linewidth=1.4, marker="o", markersize=2.5,
                label="%s (%.1f K)" % (
                    LABELS[tag],
                    statistics.fmean([
                        e["value"] for e in entries
                        if e["index"] in window_set]) if window_set else 0.0))
    frame(ax, "the near-wall cooling through the marched solution", "time (ms)",
          "$\\Delta T$ against the adiabat (K)")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.legend(fontsize=8, ncol=2)
    if stamps:
        ax.text(1000.0 * max(stamps), ax.get_ylim()[1],
                "judged window", ha="right", va="top", fontsize=8,
                color="#555555")

    fig.tight_layout()
    return save(fig, "engine-wall-bracket-effect", out_dir=out_dir)


def mode_figure(rows, out_dir):
    """the chamber's ring, and why it cannot be removed."""
    mode = gates.measure_shared_mode(rows, "mass_flow")
    series, _steps = gates.collect_mode_series(rows, "mass_flow")

    fig, axes = figure(nrows=1, ncols=2, figsize=(9.4, 4.0))
    fig.suptitle("the chamber's common mode: frequent, coherent, and not removable")

    # every leg's throat mass flow through the march, with its own mean.
    ax = axes[0]
    for tag, points in series.items():
        xs = [step for step, _ in points]
        ys = [value for _, value in points]
        ax.plot(xs, ys, marker="o", markersize=2.5, color=COLORS[tag],
                linestyle=STYLES[tag], linewidth=1.0, alpha=0.85,
                label=LABELS[tag])
    frame(ax, "the choked mass flow during the marched solution", "step",
          "mass flow (kg/s)")
    ax.legend(fontsize=8, ncol=3, loc="upper center", framealpha=0.95,
              bbox_to_anchor=(0.5, -0.16))

    # the mode's own evidence, each metric on a comparable scale.
    ax = axes[1]
    labels2 = []
    heights = []
    colours = []
    if mode.get("checked"):
        labels2.append("phase\ncoherence")
        heights.append(mode["phase_coherence"])
        colours.append(GREEN if mode.get("shares") else RED)
        labels2.append("variance\nexplained")
        heights.append(mode["adjusted_explained_fraction"])
        colours.append(PURPLE)
        labels2.append("worst\nprediction")
        heights.append(mode["weakest_cross_prediction"])
        colours.append(RED)
    if labels2:
        bars = ax.bar(range(len(labels2)), heights, color=colours, width=0.55)
        ax.set_ylim(0.0, 1.15)
        ax.axhline(1.0, color="black", linewidth=0.8, linestyle=":")
        ax.text(len(labels2) - 0.5, 1.02, "perfect", ha="right", fontsize=8)
        for rect, height in zip(bars, heights):
            ax.text(rect.get_x() + rect.get_width() / 2.0, height + 0.03,
                    "%.3f" % height, ha="center", fontsize=9)
    frame(ax, "the shared-mode diagnostic", "", "fraction")
    ax.set_xticks(range(len(labels2)))
    ax.set_xticklabels(labels2, fontsize=9)
    ax.tick_params(axis="x", pad=6)
    if mode.get("checked"):
        subtitle = (
            "period %.0f steps, %.1f cycles\n%s"
            % (mode["period_steps"], mode["cycles_covered"],
               "the frequency and phase agree"
               if mode.get("shares")
               else "but the fit cannot predict a leg"))
        ax.text(0.02, 0.97, subtitle, transform=ax.transAxes,
                ha="left", va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "#CCCCCC",
                      "boxstyle": "round,pad=0.3"})

    fig.tight_layout()
    return save(fig, "engine-wall-bracket-mode", out_dir=out_dir)


def build(runs_dir, out_dir):
    rows = gather(runs_dir)
    result = cls.assess(rows)
    return [effect_figure(rows, result, out_dir), mode_figure(rows, out_dir)]


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("runs_dir")
    parser.add_argument("out_dir")
    args = parser.parse_args(argv)
    for path in build(args.runs_dir, args.out_dir):
        print("wrote", path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

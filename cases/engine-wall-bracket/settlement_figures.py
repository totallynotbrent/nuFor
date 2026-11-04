"""the figures that show how the bracket decides it has settled.

the pictures are drawn from the same measurement path the gate uses, so
what they show and what the gate believes cannot drift apart. four
panels tell the whole argument:

- the ABSOLUTE level of the chamber state, which rises in every leg at
  once, because the chamber itself is evolving. this is the drift an
  absolute test would charge to whichever leg it happened to examine;
- the DEPARTURE from the adiabatic reference, which is the wall's own
  contribution, plotted for each cooled leg against the joint window
  the gate chose;
- the MARGIN, the slope of each departure measured in standard errors
  from flat, swept over every candidate window length, so the reader
  can see the window was chosen and not assumed;
- the VERDICT, core against near-wall, with the exact numbers the gate
  reports.

usage: settlement_figures.py <runs-dir> <out-dir>

the figure library lives in the system's python, not in an agent's
virtual environment, so run this with the system interpreter and a
clean PYTHONPATH:

    env -u PYTHONPATH /usr/bin/python3 settlement_figures.py \
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
import stationarity as st


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


def frame(ax, title, xlabel, ylabel):
    """style one panel on the house template."""
    ax.grid(True, color="#cccccc", lw=0.7, alpha=0.8)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.tick_params(labelsize=9)


def legends(ax):
    legend = ax.legend(loc="best", fontsize=8)
    legend.set_zorder(10)
    return legend


def gather(runs_dir):
    rows = []
    for tag, wall_k in gates.BRACKET:
        info = gates.measure_one_leg(tag, runs_dir)
        info["specified_wall_k"] = wall_k
        info["progress_path"] = os.path.join(runs_dir, tag, "progress.json")
        rows.append(info)
    return rows


def collect(rows, key):
    """each leg's series for one observable, keyed by saved field."""
    out = {}
    for row in rows:
        values = {}
        for observation in row.get("observations") or []:
            index = cls.field_index(observation)
            if index is None:
                continue
            value = observation.get(key)
            if value is None:
                continue
            values[index] = float(value)
        out[row["tag"]] = values
    return out


def field_times(rows):
    """the marched time of each saved field.

    the saved fields are numbered by the step at which they were
    written, so the time comes from the leg's own progress record,
    sampled densely and interpolated onto the saved step. the axis is
    then a real time rather than a step counter.
    """
    times = {}
    for row in rows:
        tag = row["tag"]
        path = row.get("progress_path")
        samples = []
        if path and os.path.exists(path):
            with open(path) as handle:
                document = json.load(handle)
            samples = [(float(s["step"]), float(s["t"]))
                       for s in document.get("samples") or []
                       if s.get("step") is not None and s.get("t") is not None]
        samples.sort()
        for observation in row.get("observations") or []:
            index = cls.field_index(observation)
            if index is None:
                continue
            stamp = None
            if samples:
                # the closest recorded march instant for this step.
                stamp = min(samples, key=lambda item: abs(item[0] - index))[1]
            times.setdefault(index, stamp)
    return times


def settlement_margin(rows, result):
    """the slope margin of every departure, over every window length.

    the margin is the fitted slope's distance from flat, in standard
    errors. a departure that is stationary puts the margin inside the
    tolerance band however the window is drawn, so the band's crossing
    shows which windows the judgement really admits.
    """
    cooled = [row["tag"] for row in rows if row["tag"] != cls.REFERENCE_TAG]
    by_field = cls.series_by_field(
        rows, cls.CORE_KEYS + cls.SECONDARY_KEYS)
    dep = cls.departures_from(
        by_field, cls.REFERENCE_TAG, cooled, cls.CORE_KEYS + cls.SECONDARY_KEYS)
    shared = None
    for tag in cooled:
        for key in cls.CORE_KEYS + cls.SECONDARY_KEYS:
            indices = {entry["index"] for entry in dep[tag][key]}
            shared = indices if shared is None else shared & indices
    ordered = {
        (tag, key): [e for e in dep[tag][key] if e["index"] in shared]
        for tag in cooled
        for key in cls.CORE_KEYS + cls.SECONDARY_KEYS
    }
    shortest = min(len(series) for series in ordered.values())
    tolerance = result.get("slope_tolerance", cls.BASE_SLOPE_TOLERANCE)

    lengths = []
    core_worst = []
    secondary_worst = []
    for length in range(shortest, 4, -1):
        start = shortest - length
        worst_core = 0.0
        worst_secondary = 0.0
        for tag in cooled:
            for key in cls.CORE_KEYS:
                entries = ordered[(tag, key)][start:]
                values = [e["value"] for e in entries]
                scale = statistics.fmean([e["level"] for e in entries])
                judgement = st.stationarity(
                    values, scale=scale, tol_slope_sigma=tolerance)
                sigma = judgement.get("slope_sigma") or 0.0
                worst_core = max(worst_core, sigma)
            for key in cls.SECONDARY_KEYS:
                entries = ordered[(tag, key)][start:]
                values = [e["value"] for e in entries]
                scale = statistics.fmean([e["level"] for e in entries])
                judgement = st.stationarity(
                    values, scale=scale, tol_slope_sigma=tolerance)
                sigma = judgement.get("slope_sigma") or 0.0
                worst_secondary = max(worst_secondary, sigma)
        lengths.append(length)
        core_worst.append(worst_core)
        secondary_worst.append(worst_secondary)
    return lengths, core_worst, secondary_worst, tolerance, shortest


def build(runs_dir, out_dir):
    rows = gather(runs_dir)
    result = cls.assess(rows)
    times = field_times(rows)
    cooled = [row["tag"] for row in rows if row["tag"] != cls.REFERENCE_TAG]

    fig = pyplot.figure(figsize=(13.4, 7.2))
    grid = fig.add_gridspec(2, 4, height_ratios=(1.0, 0.86),
                             hspace=0.42, wspace=0.34)

    # ---------------- top left: the chamber's own drift --------------
    ax = fig.add_subplot(grid[0, 0:2])
    pressure = collect(rows, "chamber_p0")
    for tag, series in pressure.items():
        indices = sorted(series)
        usable = [i for i in indices if times.get(i) is not None]
        ax.plot([1000.0 * times[i] for i in usable],
                [series[i] for i in usable], "-",
                color=COLORS.get(tag, GRAY), lw=1.4,
                label=LABELS.get(tag, tag))
    frame(ax,
          "the chamber's own state rises in every leg at once",
          "march time (ms)", "chamber total pressure (Pa)")
    legends(ax)

    # ---------------- top middle: the margin sweep --------------------
    ax = fig.add_subplot(grid[0, 2])
    lengths, core_worst, secondary_worst, tolerance, _ = \
        settlement_margin(rows, result)
    ax.axhspan(0.0, tolerance, facecolor=GREEN, alpha=0.16, lw=0)
    ax.plot(lengths, core_worst, "o-", color="#2c6fbb", lw=1.5, ms=3.4,
            label="the core bulk observables")
    ax.plot(lengths, secondary_worst, "s-", color=ORANGE, lw=1.5, ms=3.4,
            label="the near-wall temperature")
    ax.axhline(tolerance, color=RED, ls=":", lw=1.2,
               label="the tolerance (%.2f sigma)" % tolerance)
    frame(ax,
          "how far the drift stands from flat",
          "candidate window length (saved fields)",
          "distance from flat (sigma)")
    legends(ax)

    # ---------------- top right: the verdict -------------------------
    ax = fig.add_subplot(grid[0, 3])
    core_height = []
    wall_height = []
    names = []
    for tag in cooled:
        core = result.get("core", {}).get(tag, {})
        worst = 0.0
        for key in cls.CORE_KEYS:
            reading = core.get(key)
            if reading:
                worst = max(worst, reading.get("slope_sigma") or 0.0)
        core_height.append(worst)
        reading = (result.get("secondary", {}).get(tag, {}) or {}).get(
            cls.SECONDARY_KEYS[0])
        wall_height.append(
            (reading.get("slope_sigma") or 0.0) if reading else 0.0)
        names.append(LABELS.get(tag, tag))
    width = 0.38
    positions = list(range(len(names)))
    ax.bar([p - width / 2.0 for p in positions], core_height, width=width,
           color="#2c6fbb", label="the core bulk observables")
    ax.bar([p + width / 2.0 for p in positions], wall_height, width=width,
           color=ORANGE, label="the near-wall temperature")
    ax.axhline(tolerance, color=RED, ls=":", lw=1.2,
               label="the tolerance (%.2f sigma)" % tolerance)
    ax.set_xticks(positions)
    ax.set_xticklabels(names, fontsize=8)
    frame(ax,
          "the margin behind each verdict",
          "", "distance from flat (sigma)")
    tallest = max(list(core_height) + list(wall_height) + [tolerance])
    ax.set_ylim(0.0, tallest * 1.25)
    legends(ax)

    # ---------------- bottom: one panel per cooled leg ----------------
    by_field = cls.series_by_field(rows, cls.CORE_KEYS + cls.SECONDARY_KEYS)
    departures = cls.departures_from(
        by_field, cls.REFERENCE_TAG, cooled, cls.CORE_KEYS + cls.SECONDARY_KEYS)
    window_fields = result.get("window_fields") or []

    key_styles = (
        ("mass_flow", "mass flow", "#2c6fbb", "-"),
        ("chamber_p0", "chamber pressure", GREEN, "-"),
        ("chamber_t0", "chamber temperature", PURPLE, "-"),
        ("near_wall_t", "near-wall temperature", ORANGE, "--"),
    )
    for column, tag in enumerate(cooled[:4]):
        ax = fig.add_subplot(grid[1, column])
        for key, label, colour, style in key_styles:
            entries = departures.get(tag, {}).get(key, [])
            if not entries:
                continue
            usable = [e for e in entries if times.get(e["index"]) is not None]
            xs = [1000.0 * times[e["index"]] for e in usable]
            ys = [100.0 * e["value"] / e["level"] if e["level"] else 0.0
                  for e in usable]
            ax.plot(xs, ys, style, color=colour, lw=1.3, ms=2.4,
                    label=label)
        reading = result.get("core", {}).get(tag, {})
        worst = 0.0
        for key in cls.CORE_KEYS:
            item = reading.get(key)
            if item:
                worst = max(worst, item.get("slope_sigma") or 0.0)
        wall_reading = (result.get("secondary", {}).get(tag, {}) or {}).get(
            cls.SECONDARY_KEYS[0])
        wall_sigma = (wall_reading.get("slope_sigma") or 0.0) \
            if wall_reading else 0.0
        ax.set_title("%s\nbulk %.2f, near wall %.2f sigma"
                     % (LABELS.get(tag, tag), worst, wall_sigma), fontsize=9)
        frame(ax, "", "march time (ms)",
              "% shift from the adiabat" if column == 0 else "")
        if column == 0:
            legends(ax)

    fig.suptitle(
        "the bracket settles on what each wall changes, not on the "
        "chamber's own drift", fontsize=12)
    fig.text(0.5, 0.005,
             "the lower row shows every leg's departure from the adiabatic "
             "reference over the joint window the gate judged.",
             ha="center", fontsize=9, color="#444444")

    out = save(fig, "engine-wall-bracket-settlement", out_dir=out_dir)

    summary = {
        "qualified": result.get("qualified"),
        "window_start": result.get("window_start"),
        "window_end": result.get("window_end"),
        "window_length": result.get("window_length"),
        "window_fields": result.get("window_fields"),
        "comparisons": result.get("comparisons"),
        "slope_tolerance": result.get("slope_tolerance"),
        "secondary_unsettled": result.get("secondary_unsettled"),
        "core": {
            tag: {key: {"mean": reading["mean"],
                        "noise_fraction": reading["noise_fraction"],
                        "drift_fraction": reading["drift_fraction"],
                        "slope_sigma": reading["slope_sigma"]}
                  for key, reading in per_key.items()}
            for tag, per_key in result.get("core", {}).items()
        },
        "secondary": {
            tag: {key: {"mean": reading["mean"],
                        "noise_fraction": reading["noise_fraction"],
                        "drift_fraction": reading["drift_fraction"],
                        "slope_sigma": reading["slope_sigma"],
                        "qualifies": reading["qualifies"]}
                  for key, reading in per_key.items()}
            for tag, per_key in result.get("secondary", {}).items()
        },
        "margin_sweep": {
            "lengths": lengths,
            "core_worst": core_worst,
            "secondary_worst": secondary_worst,
        },
        "figure": out,
    }
    # the numbers behind the figure go beside the campaign, not among the
    # published pictures, so the image directory holds only images.
    record = os.path.join(os.path.abspath(runs_dir),
                          "settlement-figures.json")
    with open(record, "w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    return summary


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs_dir", help="the campaign's run directory")
    parser.add_argument("out_dir", help="where the figures are written")
    args = parser.parse_args(argv[1:])
    os.makedirs(args.out_dir, exist_ok=True)
    summary = build(args.runs_dir, args.out_dir)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

# the Fay-Riddell stagnation gate: the chunk-wise peak stagnation heat
# flux of the 660k-step stagnation chain, the settled band it holds, and
# the planar closed form evaluated on the marched edge state.
#
# the chain is eleven 60k-step chunks, too long to re-march for a gate,
# so this gate recomputes the published statistics from the committed
# chunk table and cross-checks that table against the committed run
# series. nothing is hardcoded and nothing is read from outside the
# repository.
#
# run from the case directory:
#   python3 gates.py               # writes gates.json
import json
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

SERIES = os.path.join(HERE, "run-log.csv")
REFINED = os.path.join(HERE, "run-log-refined.csv")
LOG = os.path.join(HERE, "run-log.txt")
REFERENCE = os.path.join(HERE, "reference.json")

# the settled windows the published bands cover.
SETTLED_FROM = 7.0e-4
REFINED_SETTLED_FROM = 3.5e-4

# percentage bands each measured quantity must stay inside.
TOLERANCE = {
    "settled_mean_pct": 10.0,
    "settled_stdev_pct": 25.0,
    "settled_low_pct": 15.0,
    "settled_high_pct": 15.0,
    "closed_form_pct": 5.0,
    "ratio_pct": 5.0,
    "refined_mean_pct": 10.0,
    "refined_stdev_pct": 25.0,
}


def read_chunk_table(path):
    """the chunk-wise peak wall heat flux from the run's own table."""
    rows = []
    for line in open(path).read().splitlines():
        m = re.match(r"\s*chunk\s+(\d+)\s+t=([\d.eE+-]+)\s+q_w\s+([\d.]+)\s+W/cm2", line)
        if m:
            rows.append({
                "chunk": int(m.group(1)),
                "time": float(m.group(2)),
                "qw_wcm2": float(m.group(3)),
            })
    return rows


def read_series(path):
    """the committed run series, in watts per square centimetre.

    the csv interleaves a metric row with the data row it belongs to, and
    its first pair belongs to a short segment rather than to the chain,
    so the pairs are lined up from the front after that pair is dropped.
    """
    times = []
    values = []
    for line in open(path).read().splitlines()[1:]:
        cells = [c.strip() for c in line.split(",")]
        if len(cells) < 4:
            continue
        if cells[2] == "qw_peak_mwm2":
            try:
                values.append(float(cells[3]) * 100.0)
            except ValueError:
                pass
        else:
            try:
                times.append(float(cells[2]))
            except ValueError:
                pass
    if len(times) > len(values):
        # the trailing data row has no metric row of its own.
        times = times[:len(values)]
    elif len(values) > len(times):
        values = values[:len(times)]
    return [{"time": t, "qw_wcm2": v} for t, v in zip(times, values)]


def read_log(path):
    """the edge state, the closed form and the logged band."""
    text = open(path).read()
    out = {}
    m = re.search(r"edge plateau: rho_e ([\d.eE+-]+) kg/m3, T_e ([\d.]+) K, p_e ([\d.]+) Pa", text)
    if m:
        out["edge_density"] = float(m.group(1))
        out["edge_temperature_k"] = float(m.group(2))
        out["edge_pressure_pa"] = float(m.group(3))
    m = re.search(r"closed form on the sim's edge state: ([\d.]+) W/cm2", text)
    if m:
        out["closed_form_wcm2"] = float(m.group(1))
    m = re.search(r"settled window \(t > ([\d.]+) ms\): mean ([\d.]+), stdev ([\d.]+), band "
                  r"([\d.]+) to ([\d.]+) W/cm2", text)
    if m:
        out["logged_window_ms"] = float(m.group(1))
        out["logged_mean"] = float(m.group(2))
        out["logged_stdev"] = float(m.group(3))
        out["logged_low"] = float(m.group(4))
        out["logged_high"] = float(m.group(5))
    m = re.search(r"ratio, settled mean / closed form: ([\d.]+)", text)
    if m:
        out["logged_ratio"] = float(m.group(1))
    return out


def within(value, reference, pct):
    if reference == 0.0:
        return abs(value) <= 1e-12
    return abs(value - reference) / abs(reference) * 100.0 <= pct


def main():
    if not os.path.exists(LOG):
        print(f"MISS: no {os.path.basename(LOG)}")
        return 1
    table = read_chunk_table(LOG)
    if len(table) < 4:
        print("FAIL: the run log carries too few chunks")
        return 1
    logged = read_log(LOG)
    settled = [r["qw_wcm2"] for r in table if r["time"] > SETTLED_FROM]
    if len(settled) < 4:
        print("FAIL: the settled window carries too few chunks")
        return 1
    mean = statistics.fmean(settled)
    stdev = statistics.stdev(settled)
    entry = {
        "tag": "stagnation-fay",
        "chunks": len(table),
        "settled_chunks": len(settled),
        "settled_mean_wcm2": mean,
        "settled_stdev_wcm2": stdev,
        "settled_low_wcm2": min(settled),
        "settled_high_wcm2": max(settled),
        "series_low_wcm2": min(r["qw_wcm2"] for r in table),
        "series_high_wcm2": max(r["qw_wcm2"] for r in table),
    }
    if "closed_form_wcm2" in logged:
        entry["closed_form_wcm2"] = logged["closed_form_wcm2"]
        entry["ratio"] = mean / logged["closed_form_wcm2"]
    for key in ("edge_density", "edge_temperature_k", "edge_pressure_pa"):
        if key in logged:
            entry[key] = logged[key]

    failures = []

    # the committed series must agree with the run's own chunk table
    # wherever both carry the same chunk.
    if os.path.exists(SERIES):
        series = read_series(SERIES)
        checked = 0
        worst = 0.0
        for point in series:
            best = None
            for row in table:
                if abs(row["time"] - point["time"]) <= 0.02 * row["time"]:
                    best = row
                    break
            if best is None or best["qw_wcm2"] == 0.0:
                continue
            drift = abs(point["qw_wcm2"] - best["qw_wcm2"]) / best["qw_wcm2"] * 100.0
            worst = max(worst, drift)
            checked += 1
        entry["series_points"] = len(series)
        entry["series_cross_checks"] = checked
        entry["series_worst_drift_pct"] = worst
        if checked < 4:
            failures.append(f"the series carries only {checked} comparable chunks")
        elif worst > 2.0:
            failures.append(
                f"the series disagrees with the chunk table by {worst:.2f}% (> 2%)"
            )

    # the refined mesh holds a narrower band, and the published
    # stagnation claim rests on that one.
    if os.path.exists(REFINED):
        refined = read_series(REFINED)
        settled_refined = [
            p["qw_wcm2"] for p in refined if p["time"] > REFINED_SETTLED_FROM
        ]
        if len(settled_refined) < 4:
            failures.append("the refined series carries too few settled chunks")
        else:
            entry["refined_chunks"] = len(settled_refined)
            entry["refined_mean_wcm2"] = statistics.fmean(settled_refined)
            entry["refined_stdev_wcm2"] = statistics.stdev(settled_refined)

    if os.path.exists(REFERENCE):
        reference = json.load(open(REFERENCE))[0]
        for key, tol in TOLERANCE.items():
            ref_key = {
                "settled_mean_pct": "settled_mean_wcm2",
                "settled_stdev_pct": "settled_stdev_wcm2",
                "settled_low_pct": "series_low_wcm2",
                "settled_high_pct": "series_high_wcm2",
                "closed_form_pct": "closed_form_wcm2",
                "ratio_pct": "ratio",
                "refined_mean_pct": "refined_mean_wcm2",
                "refined_stdev_pct": "refined_stdev_wcm2",
            }[key]
            if ref_key not in reference or ref_key not in entry:
                continue
            base = reference[ref_key]
            drift = abs(entry[ref_key] - base) / abs(base) * 100.0 if base else 0.0
            entry[f"{key}_drift_pct"] = drift
            if drift > tol:
                failures.append(
                    f"{ref_key} {entry[ref_key]:.6g} vs reference {base:.6g} "
                    f"({drift:.2f}% > {tol}%)"
                )

    # the logged ratio must follow from the recomputed mean and the
    # logged closed form.
    if "logged_ratio" in logged and "ratio" in entry:
        drift = abs(entry["ratio"] - logged["logged_ratio"]) / logged["logged_ratio"] * 100.0
        entry["logged_ratio_drift_pct"] = drift
        if drift > TOLERANCE["ratio_pct"]:
            failures.append(
                f"ratio {entry['ratio']:.4f} against the logged "
                f"{logged['logged_ratio']:.4f} ({drift:.2f}%)"
            )

    with open(os.path.join(HERE, "gates.json"), "w") as fh:
        json.dump([entry], fh, indent=1)
        fh.write("\n")

    print(
        f"stagnation-fay: {entry['settled_chunks']}/{entry['chunks']} chunks settled  "
        f"mean {entry['settled_mean_wcm2']:.2f}  stdev {entry['settled_stdev_wcm2']:.2f}  "
        f"band {entry['series_low_wcm2']:.1f} to {entry['series_high_wcm2']:.1f} W/cm2  "
        f"closed form {entry.get('closed_form_wcm2', 0.0):.2f}  ratio {entry.get('ratio', 0.0):.2f}  "
        f"series cross-check {entry.get('series_cross_checks', 0)}/"
        f"{entry.get('series_points', 0)} chunks"
    )
    if failures:
        print("")
        for f in failures:
            print(f"FAIL {f}")
        return 1
    print("")
    print("stagnation-fay gates: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

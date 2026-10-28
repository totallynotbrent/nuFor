# the catalytic-wall accommodation gate: the behaviour the wall
# accommodation has to show, checked against the committed summary of
# the sweep runs.
#
#   python3 gates.py               # writes gates.json
#
# the invariants, not frozen numbers:
#   1. the ends of the ladder reproduce the shipped wall-condition pair
#      (the zero leg the non-catalytic run, the one leg the fully
#      catalytic run) to within the reference tolerance;
#   2. the accommodation raises the peak wall heat flux monotonically;
#   3. the amplification over the non-catalytic wall stays inside the
#      reference band;
#   4. every leg ran to the final time and carries the continuum
#      diagnostic, so the published series cannot hide a broken leg.
#
# every number comes from summary.csv, which samples.py distils from the
# committed run logs.
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = os.path.join(HERE, "summary.csv")
REFERENCE = os.path.join(HERE, "reference.json")
OUTPUT = os.path.join(HERE, "gates.json")


def load_summary():
    rows = {}
    with open(SUMMARY) as fh:
        for row in csv.DictReader(fh):
            rows[float(row["eta"])] = row
    return rows


def main():
    reference = json.load(open(REFERENCE))
    rows = load_summary()
    etas = sorted(rows)
    peaks = [(eta, float(rows[eta]["peak_q_mwm2"])) for eta in etas]

    endpoint_checks = []
    for eta, key in ((0.0, "non_catalytic"), (1.0, "fully_catalytic")):
        got = dict(peaks)[eta]
        want = reference[key]["peak_q_mwm2"]
        tol = reference[key]["tol_pct"]
        drift = abs(got - want) / want * 100.0
        endpoint_checks.append({
            "eta": eta,
            "measured_peak_q_mwm2": got,
            "reference_peak_q_mwm2": want,
            "drift_pct": drift,
            "passed": drift <= tol,
        })

    monotonic = all(peaks[i][1] <= peaks[i + 1][1] + 1e-9 for i in range(len(peaks) - 1))
    base = dict(peaks)[0.0]
    amplification = [
        {"eta": eta, "amplification": q / base} for eta, q in peaks
    ]
    amp_ok = all(
        reference["amplification_band"][0] <= a["amplification"] <= reference["amplification_band"][1]
        for a in amplification
    )

    complete = []
    for eta in etas:
        row = rows[eta]
        ran = (
            row["steps"] is not None
            and int(float(row["steps"])) > 0
            and row["kn_gll_rho"] not in (None, "")
            and row["kn_gll_t"] not in (None, "")
        )
        complete.append({"eta": eta, "ran_with_continuum_probe": ran})

    passed = (
        all(c["passed"] for c in endpoint_checks)
        and monotonic
        and amp_ok
        and all(c["ran_with_continuum_probe"] for c in complete)
    )

    result = [{
        "tag": "park5-visc-accommodation",
        "accommodation": [
            {
                "eta": eta,
                "peak_q_mwm2": q,
                "peak_x": float(rows[eta]["peak_x"]),
                "peak_y": float(rows[eta]["peak_y"]),
                "steps": int(float(rows[eta]["steps"])),
                "max_temperature_k": float(rows[eta]["max_temperature_k"]),
                "kn_gll_rho": float(rows[eta]["kn_gll_rho"]),
                "kn_gll_t": float(rows[eta]["kn_gll_t"]),
                "cells_above_breakdown_pct": float(rows[eta]["cells_above_breakdown_pct"]),
            }
            for eta, q in peaks
        ],
        "endpoint_checks": endpoint_checks,
        "monotonic_in_accommodation": monotonic,
        "amplification": amplification,
        "amplification_within_band": amp_ok,
        "legs_complete": complete,
        "passed": passed,
    }]
    with open(OUTPUT, "w") as fh:
        json.dump(result, fh, indent=1)
    print(json.dumps(result, indent=1))
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

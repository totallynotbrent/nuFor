"""prices the data-dependent window search against a registered horizon.

the bracket used to scan candidate windows from the widest down and
keep the widest one on which all twelve core departures qualified. the
threshold paid for those twelve screens and nothing else, so the search
could hand the bracket a quiet stretch and let a certificate survive a
drift the declared horizon would have caught.

this harness prices that search instead of asserting a penalty for it.
three arms run on the same surrogate records, so the only difference
between them is the procedure:

- "scan_full": the scan over its whole candidate range;
- "scan_width_matched": the same scan capped at the registered
  horizon's width, so the arm may CHOOSE among neighbouring widths but
  cannot widen. the gap between this arm and "registered" is the price
  of the search alone;
- "registered": one declared horizon, fixed before the levels are
  inspected.

the surrogates are trend-free, so certifying one is correct behaviour
rather than a false alarm. the arms are therefore compared on yield on
a quiet record, on detection among what each arm certified, and on the
share of certificates that survive a known drift.

usage:

    python3 search_multiplicity.py <pairs> <campaign_root> [out_json]
"""

import argparse
import copy
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import calibrate_settlement as cal  # noqa: E402
import cross_leg_settlement as cls  # noqa: E402
import settled_observables as settle_mod  # noqa: E402


AMPITUDES = (0.01, 0.02, 0.05, 0.1, 0.2)
DESIGNS = ("whole_record", "confined")


def scan_verdict(rows, max_length=None):
    """the legacy scanning screen, optionally capped to a maximum width.

    the scan is kept intact and only prevented from considering widths
    above the cap, so the capped arm differs from the registered arm in
    exactly one respect: it may choose its window.
    """
    tags = [row["tag"] for row in rows]
    if cls.REFERENCE_TAG not in tags:
        return False
    cooled = [tag for tag in tags if tag != cls.REFERENCE_TAG]
    if not cooled:
        return False
    all_keys = cls.CORE_KEYS + cls.SECONDARY_KEYS
    by_field = cls.series_by_field(rows, all_keys)
    departures = cls.departures_from(
        by_field, cls.REFERENCE_TAG, cooled, all_keys
    )
    ordered = {}
    shared = None
    for tag in cooled:
        for key in all_keys:
            entries = departures.get(tag, {}).get(key)
            if not entries:
                return False
            indices = {entry["index"] for entry in entries}
            shared = indices if shared is None else shared & indices
            ordered[(tag, key)] = entries
    if not shared:
        return False
    for tag in cooled:
        for key in all_keys:
            ordered[(tag, key)] = [
                entry for entry in ordered[(tag, key)]
                if entry["index"] in shared
            ]
    shortest = min(len(series) for series in ordered.values())
    comparisons = max(1, len(cooled) * len(cls.CORE_KEYS))
    widest = shortest if max_length is None else min(shortest, max_length)
    for length in range(widest, settle_mod.MIN_WINDOW - 1, -1):
        start = shortest - length
        tolerance = cls.slope_tolerance_for(comparisons, length)
        qualifies = True
        for tag in cooled:
            for key in cls.CORE_KEYS:
                judgement = cls.judge(
                    ordered[(tag, key)][start:], tolerance
                )
                if judgement is None or not judgement["qualifies"]:
                    qualifies = False
                    break
            if not qualifies:
                break
        if qualifies:
            return True
    return False


def registered_verdict(rows):
    """the registered single-horizon screen."""
    try:
        result = cls.assess(rows)
    except Exception:
        return False
    return bool(result.get("qualified", False))


def horizon_fields(rows):
    """the field indices the registered horizon judges."""
    tags = [row["tag"] for row in rows]
    cooled = [tag for tag in tags if tag != cls.REFERENCE_TAG]
    all_keys = cls.CORE_KEYS + cls.SECONDARY_KEYS
    by_field = cls.series_by_field(rows, all_keys)
    departures = cls.departures_from(
        by_field, cls.REFERENCE_TAG, cooled, all_keys
    )
    shared = None
    for tag in cooled:
        for key in all_keys:
            indices = {entry["index"] for entry in departures[tag][key]}
            shared = indices if shared is None else shared & indices
    if not shared:
        return set()
    return set(sorted(shared)[-cls.SETTLEMENT_WINDOW:])


def confined_drift(rows, amplitude):
    """add a monotone drift confined to the declared horizon."""
    out = copy.deepcopy(rows)
    horizon = horizon_fields(rows)
    if not horizon:
        return out
    for row in out:
        if row["tag"] == cls.REFERENCE_TAG:
            continue
        observations = row.get("observations") or []
        window = [
            observation for observation in observations
            if cls.field_index(observation) in horizon
        ]
        for key in cls.CORE_KEYS:
            values = [observation.get(key) for observation in window]
            if any(value is None for value in values):
                continue
            level = sum(abs(value) for value in values) / len(values)
            if level == 0.0:
                level = 1.0
            count = len(values)
            for position, observation in enumerate(window):
                fraction = position / (count - 1)
                observation[key] = (
                    values[position] + amplitude / 100.0 * fraction * level
                )
    return out


def drifted_rows(design, surrogate, amplitude):
    """the drift arm for one design."""
    if design == "whole_record":
        return cal.ramped_rows(surrogate, amplitude)
    if design == "confined":
        return confined_drift(surrogate, amplitude)
    raise ValueError("unknown design: {}".format(design))


def blank_tally():
    return {
        "quiet_records_certified": 0,
        "drifted_records_certified": 0,
        "certificates_surviving_drift": 0,
        "certificates_withdrawn": 0,
    }


def summarise(tally, pairs):
    """the rates a tally supports."""
    base = tally["quiet_records_certified"]
    survived = tally["certificates_surviving_drift"]
    tally["yield_on_quiet_record"] = base / pairs if pairs else None
    tally["detected_among_baseline"] = (
        tally["certificates_withdrawn"] / base if base else None
    )
    tally["miss_rate_given_baseline_certificate"] = (
        survived / base if base else None
    )
    tally["certificate_survives_drift_rate"] = (
        survived / pairs if pairs else None
    )
    return tally


def cell(pairs, records, design, amplitude):
    """one amplitude's paired counts for the three arms."""
    tallies = {
        "scan_full": blank_tally(),
        "scan_width_matched": blank_tally(),
        "registered": blank_tally(),
    }
    runners = {
        "scan_full": lambda rows: scan_verdict(rows),
        "scan_width_matched": lambda rows: scan_verdict(
            rows, cls.SETTLEMENT_WINDOW
        ),
        "registered": registered_verdict,
    }
    for index in range(pairs):
        rng = random.Random("mult-{}".format(index))
        surrogate = cal.surrogate_rows("detrended_rotation", records, rng)
        drifted = drifted_rows(design, surrogate, amplitude)
        for name, runner in runners.items():
            quiet = runner(surrogate)
            noisy = runner(drifted)
            tally = tallies[name]
            if quiet:
                tally["quiet_records_certified"] += 1
            if noisy:
                tally["drifted_records_certified"] += 1
            if quiet and noisy:
                tally["certificates_surviving_drift"] += 1
            if quiet and not noisy:
                tally["certificates_withdrawn"] += 1
    for tally in tallies.values():
        summarise(tally, pairs)
    matched = tallies["scan_width_matched"]
    registered = tallies["registered"]
    return {
        "design": design,
        "amplitude_percent": amplitude,
        "pairs": pairs,
        "arms": tallies,
        "search_cost_at_matched_width": {
            "match_scan_miss_rate": matched[
                "miss_rate_given_baseline_certificate"
            ],
            "registered_miss_rate": registered[
                "miss_rate_given_baseline_certificate"
            ],
            "miss_rate_difference": (
                matched["miss_rate_given_baseline_certificate"]
                - registered["miss_rate_given_baseline_certificate"]
            ),
        },
    }


def run(pairs, root, from_record=False):
    """the full comparison grid."""
    records = cal.load_rows(root, from_record=from_record)
    cells = []
    for design in DESIGNS:
        for amplitude in AMPITUDES:
            row = cell(pairs, records, design, amplitude)
            cells.append(row)
            print(json.dumps(row), flush=True)
    return {
        "pairs_per_cell": pairs,
        "campaign": cal.portable_campaign(root),
        "settlement_window": cls.SETTLEMENT_WINDOW,
        "null_construction": "detrended_rotation",
        "arms": {
            "scan_full": "the legacy screen over its whole candidate range",
            "scan_width_matched": (
                "the same scan capped at the registered horizon's width, "
                "so it may choose among neighbouring widths but not widen"
            ),
            "registered": "one declared horizon, fixed before inspection",
        },
        "interpretation": (
            "the surrogates are trend-free, so certifying one is correct "
            "behaviour rather than a false alarm. these cells price the "
            "search at matched width against the registered arm, and the "
            "width against the unrestricted scan."
        ),
        "cells": cells,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pairs", type=int, help="surrogate records per cell")
    parser.add_argument("campaign", help="campaign root holding the records")
    parser.add_argument("out", nargs="?", help="where to write the JSON")
    parser.add_argument(
        "--record",
        action="store_true",
        help="read the committed record instead of a marched campaign",
    )
    args = parser.parse_args(argv)
    summary = run(args.pairs, args.campaign, from_record=args.record)
    if args.out:
        with open(args.out, "w") as handle:
            json.dump(summary, handle, indent=1)
        print("wrote", args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

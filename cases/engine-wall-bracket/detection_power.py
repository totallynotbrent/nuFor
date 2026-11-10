"""dedicated detection-power study for the bracket settlement screen.

each amplitude gets its own random stream, and every arm is a paired
comparison of one surrogate against the same surrogate with a drift
added. the SURROGATE is drawn from the trend-free null, so the only
physical difference between the arms is the injected drift.

two drift designs are measured:

- "whole_record" spreads the drift across the complete record;
- "confined" confines it to the declared settlement horizon.

the metrics are kept separate on purpose. "detected_among_baseline"
measures the screen's sensitivity among pairs that certified before the
drift arrived. "certificate_survives_drift_rate" measures the
unconditional probability that a certificate survives, including pairs
that did not certify to begin with. the former answers whether the
screen notices; the latter answers how often a bad bracket escapes.

usage:

    python3 detection_power.py <pairs> <campaign_root> [out_json] \
        [--null=detrended_rotation|raw_shift] \
        [--design=whole_record|confined|both]
"""

import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import calibrate_settlement as cal  # noqa: E402
import cross_leg_settlement as cls  # noqa: E402


AMPS = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)


def wilson(count, total, z=1.959964):
    """wilson score interval for a proportion."""
    if total == 0:
        return None, None
    phat = count / total
    denom = 1.0 + z * z / total
    centre = (phat + z * z / (2 * total)) / denom
    half = (
        z
        * ((phat * (1 - phat) / total + z * z / (4 * total * total)) ** 0.5)
        / denom
    )
    return centre - half, centre + half


def shared_tail_indices(rows):
    """the ordered fields covered by the declared settlement horizon."""
    tags = [row["tag"] for row in rows]
    cooled = [tag for tag in tags if tag != cls.REFERENCE_TAG]
    all_keys = cls.CORE_KEYS + cls.SECONDARY_KEYS
    by_field = cls.series_by_field(rows, all_keys)
    dep = cls.departures_from(by_field, cls.REFERENCE_TAG, cooled, all_keys)
    shared = None
    for tag in cooled:
        for key in all_keys:
            indices = {entry["index"] for entry in dep[tag][key]}
            shared = indices if shared is None else shared & indices
    if not shared:
        return []
    return sorted(shared)[-cls.SETTLEMENT_WINDOW:]


def confined_drift_rows(rows, amplitude, horizon_indices):
    """add a monotone drift across only the settlement horizon."""
    import copy

    out = copy.deepcopy(rows)
    inside = set(horizon_indices)
    for row in out:
        if row["tag"] == cls.REFERENCE_TAG:
            continue
        observations = row.get("observations") or []
        window = [
            observation for observation in observations
            if cls.field_index(observation) in inside
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
                    values[position]
                    + amplitude / 100.0 * fraction * level
                )
    return out


def drift_rows(design, surrogate, amplitude, horizon_indices):
    """the drift arm for one design."""
    if design == "whole_record":
        return cal.ramped_rows(surrogate, amplitude)
    if design == "confined":
        return confined_drift_rows(
            surrogate, amplitude, horizon_indices
        )
    raise ValueError("unknown drift design: {}".format(design))


def study(pairs, root, null_kind, design, from_record=False):
    """the paired power study for one null construction and design."""
    rows = cal.load_rows(root, from_record=from_record)
    horizon_indices = shared_tail_indices(rows)
    results = []
    for amplitude in AMPS:
        baseline = 0
        survived = 0
        withdrew = 0
        appeared = 0
        neither = 0
        for index in range(pairs):
            rng = random.Random(
                "{}-{}-{}-{}".format(
                    design, null_kind, amplitude, index
                )
            )
            surrogate = cal.surrogate_rows(null_kind, rows, rng)
            q_none = cal.assessed(surrogate)[0]
            q_drift = cal.assessed(
                drift_rows(design, surrogate, amplitude, horizon_indices)
            )[0]
            if q_none:
                baseline += 1
            if q_drift:
                survived += 1
            if q_none and not q_drift:
                withdrew += 1
            elif q_drift and not q_none:
                appeared += 1
            elif not q_none and not q_drift:
                neither += 1
        conditional_lo, conditional_hi = (
            wilson(withdrew, baseline) if baseline else (None, None)
        )
        survival_lo, survival_hi = wilson(survived, pairs)
        results.append({
            "design": design,
            "null_construction": null_kind,
            "amplitude_percent": amplitude,
            "pairs": pairs,
            "baseline_certificates": baseline,
            "drift_arm_certificates": survived,
            "certificate_survives_drift_rate": survived / pairs if pairs else None,
            "certificate_survives_drift_interval": [survival_lo, survival_hi],
            "detected_among_baseline": withdrew / baseline if baseline else None,
            "detected_among_baseline_interval": [
                conditional_lo, conditional_hi
            ],
            "certificate_withdrawn": withdrew,
            "certificate_appeared": appeared,
            "neither_arm_certified": neither,
        })
        print(json.dumps(results[-1]), flush=True)
    return results


def main(argv):
    args = [argument for argument in argv[1:] if not argument.startswith("--")]
    flags = [argument for argument in argv[1:] if argument.startswith("--")]
    root = args[1] if len(args) > 1 else None
    out = args[2] if len(args) > 2 else None
    pairs = int(args[0]) if args else 800
    null_kind = "detrended_rotation"
    design = "both"
    from_record = False
    for flag in flags:
        if flag.startswith("--null="):
            null_kind = flag.split("=", 1)[1]
        elif flag.startswith("--design="):
            design = flag.split("=", 1)[1]
        elif flag == "--record":
            from_record = True
    if root is None:
        print(
            "usage: detection_power.py <pairs> <campaign_root> [out_json]"
            " [--null=detrended_rotation|raw_shift]"
            " [--design=whole_record|confined|both]"
        )
        return 2
    designs = ("whole_record", "confined") if design == "both" else (design,)
    results = []
    for selected in designs:
        results.extend(study(pairs, root, null_kind, selected, from_record))
    if out:
        with open(out, "w") as handle:
            json.dump({
                "pairs_per_amplitude": pairs,
                "campaign": cal.portable_campaign(root),
                "settlement_window": cls.SETTLEMENT_WINDOW,
                "arms": results,
            }, handle, indent=1)
        print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))

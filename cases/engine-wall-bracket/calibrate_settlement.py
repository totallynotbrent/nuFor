"""calibrate the bracket's settlement screen against a surrogate null.

what this measures

the bracket judges a DECLARED final horizon on the twelve core
departures and then reports that the bracket settled. the horizon does
not move with the levels, so the verdict cannot benefit from choosing a
favourable stretch after seeing the data.

there are two ways for that procedure to be wrong, and they were
originally confused for one another.

the MISS is the dangerous one: a real drift in the wall's effect slips
through and the bracket is declared settled. it is measured by the
paired detection study below, which assesses each surrogate twice, once
as it stands and once with a known ramp added, so the ONLY difference
between the arms is the ramp.

the FALSE ALARM is the awkward one: a bracket whose wall effect really
is flat is thrown out because the search happened to dislike some
stretch of its record. on a FLAT null, certifying is the CORRECT
outcome, so the rate of interest here is the complement of the
certification rate. a screen that throws out flat brackets rarely is
not conservative, it is unusable.

the surrogate constructions

"raw_shift" rotates each cooled leg's own values. it does NOT produce a
valid null for the screen's question, because rotating a drifting series
winds the drift onto a seam instead of removing it: the surrogate still
carries the trend the null is meant to exclude. it is kept so the older
numbers can be reproduced for comparison, and it must never be quoted as
the screen's error rate.

"detrended_rotation" takes each leg's DEPARTURE from the adiabatic
reference, fits and removes its linear trend, rotates the residual, and
rebuilds the leg so its departure is the mean departure plus the rotated
residual. the trend is gone by construction, the neighbour-to-neighbour
structure of the residual survives, and the reference keeps the common
trajectory the subtraction cancels.

the construction changes the answer by an order of magnitude, so no
single surrogate number can be quoted as the error rate without naming
the construction that produced it.

how to read the output

"flat_null_false_alarm_rate" is the probability the procedure throws
out a bracket whose wall effect really is flat, measured on the flat
null. "flat_null_certified_frequency" is its complement, and on a
properly flat null it SHOULD be high: certifying a flat bracket is the
desired behaviour.

"detection" holds the paired study's per-amplitude results. a pair
where the certificate is withdrawn when the ramp is added is evidence
the screen heard the ramp. the smallest ramp whose withdrawal rate
clears the stated target is the finest drift the screen resolves, which
bounds the certificate in the honest direction.
"""

import argparse
import copy
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import cross_leg_settlement as cls


# the number of surrogate records the calibration builds.
DEFAULT_TRIALS = 400

# the seed the calibration starts from, so a run repeats exactly.
DEFAULT_SEED = 20261102

# the confidence the detection study aims for.
DETECTION_TARGET = 0.80

# the ramps the detection study injects, as a percentage of the
# observable's own level.
DETECTION_AMPLICADES = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2)

# how many baseline certifications an amplitude arm needs before its
# conditional detection rate means anything at all.
MINIMUM_DETECTION_BASELINE = 15


def portable_campaign(root):
    """the campaign path as it should be written into a shipped artefact.

    an absolute path would tie the artefact to one machine and name the
    owner's home directory, so the record folder is named relatively.
    """
    text = os.path.abspath(str(root)).replace(os.sep, "/")
    if "/record" in text:
        return text.rsplit("/", 1)[-1]
    return os.path.basename(text.rstrip("/")) or text


def load_rows(campaign_root, from_record=False):
    """the analyser's own rows for a campaign, loaded through the gate.

    the gate reads the saved fields and measures the observables, so the
    calibration sees exactly the numbers the published screen saw. the
    record path reads the same measurements from the committed record,
    so the studies can be reproduced from what the repository ships.
    """
    import gates

    rows = []
    for tag, wall_k in gates.bracket.BRACKET:
        if from_record:
            info = gates.measure_one_leg_from_record(tag)
        else:
            info = gates.measure_one_leg(tag, campaign_root)
        info["specified_wall_k"] = wall_k
        rows.append(info)
    return [row for row in rows if row.get("tag")]



def wilson_interval(count, total, z=1.959964):
    """wilson score interval for a proportion."""
    if not total:
        return None, None
    phat = count / total
    denominator = 1.0 + z * z / total
    centre = (phat + z * z / (2 * total)) / denominator
    half = (
        z
        * ((phat * (1 - phat) / total + z * z / (4 * total * total)) ** 0.5)
        / denominator
    )
    return centre - half, centre + half

def shift_cycle(values, amount):
    """a circular shift of a series, preserving its cyclic order."""
    if not values:
        return list(values)
    amount = amount % len(values)
    if amount == 0:
        amount = 1
    rotated = list(values)
    for _ in range(amount):
        rotated = rotated[-1:] + rotated[:-1]
    return rotated


def shifted_rows(rows, rng):
    """the rows with each cooled leg's values circularly shifted.

    the reference leg is untouched, because its common trajectory is
    what the departure subtracts. scrubbing it too would destroy the
    cancellation the judgement depends on.
    """
    out = copy.deepcopy(rows)
    for row in out:
        if row["tag"] == cls.REFERENCE_TAG:
            continue
        observations = row.get("observations") or []
        for key in cls.CORE_KEYS:
            values = [o.get(key) for o in observations]
            if any(value is None for value in values):
                continue
            amount = rng.randrange(1, max(len(values), 2))
            rolled = shift_cycle(values, amount)
            for observation, value in zip(observations, rolled):
                observation[key] = value
    return out


def detrended_departure_rows(rows, rng):
    """surrogate records whose wall departures carry no secular trend.

    the judgement is made on each leg's DEPARTURE from the adiabatic
    reference, so the null of interest is "the departure has no trend".
    a circular shift of the raw leg values does not produce that null:
    shifting a drifting series keeps the drift in the record, just
    wound onto a seam, so the surrogate still carries the trend the
    null is meant to exclude. shifting therefore measures the screen's
    reaction to a wandering series rather than to a stationary one.

    this constructs the null in the space the screen actually uses. for
    each cooled leg and core observable it takes the departure, fits and
    removes its linear trend, rotates the remaining residuals along
    their cyclic order, and rebuilds the leg so that its departure from
    the untouched reference is the mean departure plus the rotated
    residuals. the trend is gone by construction, the residual's
    neighbour-to-neighbour structure survives the rotation, and the
    reference keeps the common trajectory the subtraction cancels.
    """
    import statistics as stats

    out = copy.deepcopy(rows)

    reference_series = {}
    for row in rows:
        if row["tag"] != cls.REFERENCE_TAG:
            continue
        for observation in row.get("observations") or []:
            index = cls.field_index(observation)
            if index is None:
                continue
            for key in cls.CORE_KEYS:
                value = observation.get(key)
                if value is not None:
                    reference_series.setdefault(key, {})[index] = float(value)

    for row in out:
        if row["tag"] == cls.REFERENCE_TAG:
            continue
        observations = row.get("observations") or []
        for key in cls.CORE_KEYS:
            reference = reference_series.get(key) or {}
            indexed = []
            for observation in observations:
                index = cls.field_index(observation)
                value = observation.get(key)
                if index is None or value is None or index not in reference:
                    continue
                indexed.append((observation, index, float(value)))
            if len(indexed) < 4:
                continue
            departures = [
                value - reference[index] for _observation, index, value in indexed
            ]
            count = len(departures)
            mean_x = (count - 1) / 2.0
            mean_y = stats.fmean(departures)
            sxx = sum((i - mean_x) ** 2 for i in range(count))
            if sxx <= 0.0:
                continue
            sxy = sum(
                (i - mean_x) * (departures[i] - mean_y) for i in range(count)
            )
            slope = sxy / sxx
            intercept = mean_y - slope * mean_x
            residuals = [
                departures[i] - (intercept + slope * i) for i in range(count)
            ]
            amount = rng.randrange(1, max(count, 2))
            rotated = shift_cycle(residuals, amount)
            for position, (observation, index, _value) in enumerate(indexed):
                observation[key] = reference[index] + mean_y + rotated[position]
    return out


def ramped_rows(rows, amplitude):
    """the rows with a known ramp added to each cooled leg.

    the ramp grows linearly across the record and is proportional to
    the observable's own level, so the injected size means the same
    thing in every leg and every observable.
    """
    out = copy.deepcopy(rows)
    for row in out:
        if row["tag"] == cls.REFERENCE_TAG:
            continue
        observations = row.get("observations") or []
        count = len(observations)
        if count < 2:
            continue
        for key in cls.CORE_KEYS:
            values = [o.get(key) for o in observations]
            if any(value is None for value in values):
                continue
            level = sum(abs(value) for value in values) / count
            if level == 0.0:
                level = 1.0
            for position, observation in enumerate(observations):
                fraction = position / (count - 1)
                observation[key] = values[position] + (
                    amplitude / 100.0 * fraction * level
                )
    return out


def assessed(rows):
    """the procedure's verdict on the rows, and the window it chose."""
    result = cls.assess(rows)
    return (
        bool(result.get("qualified")),
        result.get("window_length"),
        result.get("slope_tolerance"),
    )


# the surrogate constructions the calibration may draw its null from.
#
# "raw_shift" rotates each cooled leg's own values. it is NOT a valid
# null for the screen's own question, because rotating a drifting series
# winds the drift onto a seam instead of removing it: the surrogate
# still carries a trend the null is meant to exclude. it is kept only
# so the older number can be reproduced for comparison.
#
# "detrended_rotation" removes each departure's linear trend first and
# rotates the residual, which is the null the screen's question really
# asks, though how much structure that surrogate should retain is still
# under study.
NULL_CONSTRUCTIONS = ("raw_shift", "detrended_rotation")


def surrogate_rows(kind, rows, rng):
    """the surrogate record for one null construction."""
    if kind == "raw_shift":
        return shifted_rows(rows, rng)
    if kind == "detrended_rotation":
        return detrended_departure_rows(rows, rng)
    raise ValueError("unknown null construction: {}".format(kind))


def calibrate(campaign_root, trials=DEFAULT_TRIALS, seed=DEFAULT_SEED,
              null_kind="detrended_rotation", from_record=False):
    """the full calibration of the settlement screen."""
    rows = load_rows(campaign_root, from_record=from_record)
    rng = random.Random(seed)

    observed_qualified, observed_length, observed_tolerance = assessed(rows)

    null_qualified = 0
    null_lengths = []
    null_tolerances = []
    for _ in range(trials):
        qualified, length, tolerance = assessed(
            surrogate_rows(null_kind, rows, rng)
        )
        if qualified:
            null_qualified += 1
            null_lengths.append(length)
            if tolerance is not None:
                null_tolerances.append(tolerance)

    detection = []
    minimum_detectable = None
    detection_trials = max(trials // 4, 40)
    for amplitude in DETECTION_AMPLICADES:
        # the two arms share every surrogate, so the only difference
        # between them is the injected ramp.
        certified_without = 0
        certified_with = 0
        withdrew = 0
        appeared = 0
        lengths = []
        for index in range(detection_trials):
            trial_rng = random.Random(
                "{}-{}-{}".format(DEFAULT_SEED, amplitude, index)
            )
            surrogate = surrogate_rows(null_kind, rows, trial_rng)
            qualified_without, _length, _tol = assessed(surrogate)
            prepared = ramped_rows(surrogate, amplitude)
            qualified_with, length, _tol2 = assessed(prepared)
            if qualified_without:
                certified_without += 1
            if qualified_with:
                certified_with += 1
                lengths.append(length)
            if qualified_without and not qualified_with:
                withdrew += 1
            elif qualified_with and not qualified_without:
                appeared += 1
        without_rate = certified_without / detection_trials if detection_trials else 0.0
        with_rate = certified_with / detection_trials if detection_trials else 0.0
        discordant = withdrew + appeared
        # the informative denominator is the surrogates the screen
        # would otherwise certify. of those, the ones that lose the
        # certificate when the ramp arrives are the detections.
        conditional = withdrew / certified_without if certified_without else None
        detection.append({
            "amplitude_percent": amplitude,
            "trials": detection_trials,
            "certified_without_ramp": certified_without,
            "certified_without_ramp_rate": without_rate,
            "certified_with_ramp": certified_with,
            "certified_with_ramp_rate": with_rate,
            "certificate_withdrawn": withdrew,
            "certificate_appeared": appeared,
            "discordant_pairs": discordant,
            "conditional_detection_rate": conditional,
            "chosen_lengths": sorted(set(lengths)),
        })
        usable = (
            conditional is not None
            and certified_without >= MINIMUM_DETECTION_BASELINE
            and conditional >= DETECTION_TARGET
        )
        if usable and minimum_detectable is None:
            minimum_detectable = amplitude

    return {
        "campaign": portable_campaign(campaign_root),
        "seed": seed,
        "trials": trials,
        "null_construction": null_kind,
        "observed_qualified": observed_qualified,
        "observed_window_length": observed_length,
        "observed_slope_tolerance": observed_tolerance,
        "family_alpha": cls.FAMILY_ALPHA,
        "flat_null_certified_count": null_qualified,
        "flat_null_certified_frequency": (
            null_qualified / trials if trials else None
        ),
        "flat_null_false_alarm_rate": (
            (trials - null_qualified) / trials if trials else None
        ),
        "flat_null_false_alarm_interval": wilson_interval(
            trials - null_qualified, trials
        ),
        "flat_null_distinct_lengths": sorted(set(null_lengths)),
        "null_tolerance_range": [
            min(null_tolerances) if null_tolerances else None,
            max(null_tolerances) if null_tolerances else None,
        ],
        "detection": detection,
        "detection_target": DETECTION_TARGET,
        "minimum_detectable_percent": minimum_detectable,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign")
    parser.add_argument("--trials", type=int, default=DEFAULT_TRIALS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--json", default=None)
    parser.add_argument(
        "--record",
        action="store_true",
        help="read the committed record instead of a marched campaign",
    )
    parser.add_argument(
        "--null",
        choices=NULL_CONSTRUCTIONS,
        default="detrended_rotation",
        help="which surrogate null to draw the false-certificate rate from",
    )
    args = parser.parse_args(argv)

    summary = calibrate(
        args.campaign,
        args.trials,
        args.seed,
        null_kind=args.null,
        from_record=args.record,
    )
    print(json.dumps(summary, indent=1))
    if args.json:
        with open(args.json, "w") as handle:
            json.dump(summary, handle, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

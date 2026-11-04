"""the bracket's settlement judgement, made on what each leg adds.

the five legs march the same chamber from the same initial condition, so
every leg carries the chamber's own evolution: the chamber warms, its
pressures drift, its mass flow walks. that evolution is COMMON. a
criterion applied to one leg's ABSOLUTE level charges the leg for a
transient it does not own, and it charges the adiabatic reference
hardest of all, because the reference draws no heat and so has the
smallest excursion to blame. that asymmetry is not a measurement, it is
an artefact of asking the wrong question, and it used to fail a bracket
whose contrasts were demonstrably steady.

the JUDGEMENT is therefore made on each leg's WALL EFFECT: its value at
a saved field LESS the adiabatic reference's value at the same field.
the reference supplies the common trajectory, the subtraction cancels
it, and what is left is the wall's contribution, which is the quantity
the bracket compares.

three things make the judgement defensible:

- the departure is judged against the LEVEL it came from, the
  reference's own magnitude, so a small difference is not flattered by
  being small;
- the window is JOINT, so a level and a difference are never quoted
  from different stretches of the march;
- the SCREEN is multiple-comparison AWARE. twelve departures are
  examined at once (four legs against the reference, three core
  observables each), so a naive two-sigma slope test would fail a
  perfectly settled bracket by chance alone. the tolerance is widened
  for that many comparisons. the widening treats the twelve screens as
  the whole family: the judgement horizon is declared in advance, so no
  window search joins the family. the residuals are serially dependent
  beyond the lag-one correction, so the achieved rejection rate is
  MEASURED on surrogate records (calibrate_settlement.py), not
  promised by the algebra.

the near-wall temperature is judged SEPARATELY from the core bulk
observables. it is a wall's sharpest signature, but it is also the
quantity that keeps moving while the bulk settles, so folding it into
the core test would hide one behind the other. the result states which
is which, and the near-wall verdict is reported even when it fails.
"""
import math
import os
import statistics

import significance as sig

import settled_observables as settle_mod
import stationarity as st


CORE_KEYS = ("mass_flow", "chamber_p0", "chamber_t0")
SECONDARY_KEYS = ("near_wall_t",)
REFERENCE_TAG = "wall-adiabatic"

# the base slope tolerance of the stationarity test, before the screen's
# multiple comparisons are accounted for.
BASE_SLOPE_TOLERANCE = 2.0

# the error budget the whole screen spends, split across the comparisons.
FAMILY_ALPHA = 0.05

# the final horizon judged on every run. the horizon is fixed before the
# levels are inspected, so the verdict cannot benefit from searching the
# record for a favourable suffix.
SETTLEMENT_WINDOW = 14


def field_index(observation):
    """the saved field's ordinal, taken from its file name."""
    name = os.path.basename(observation.get("path", ""))
    try:
        return int(name.split("-")[1].split(".")[0])
    except (IndexError, ValueError):
        return None


def series_by_field(rows, keys):
    """each observable's series per leg, keyed by saved field."""
    out = {}
    for row in rows:
        tag = row["tag"]
        for observation in row.get("observations") or []:
            index = field_index(observation)
            if index is None:
                continue
            for key in keys:
                value = observation.get(key)
                if value is None:
                    continue
                out.setdefault(key, {}).setdefault(tag, {})[index] = float(value)
    return out


def departures_from(by_field, reference, tags, keys):
    """each leg's departure from the reference, over the shared fields."""
    out = {}
    for tag in tags:
        per_key = {}
        for key in keys:
            reference_values = (by_field.get(key) or {}).get(reference) or {}
            leg_values = (by_field.get(key) or {}).get(tag) or {}
            shared = sorted(set(reference_values) & set(leg_values))
            if not shared:
                continue
            per_key[key] = [
                {
                    "index": index,
                    "value": leg_values[index] - reference_values[index],
                    "level": abs(reference_values[index]),
                }
                for index in shared
            ]
        out[tag] = per_key
    return out


def judge(entries, slope_tolerance=BASE_SLOPE_TOLERANCE):
    """one departure over one window, judged against the level it came from."""
    if not entries:
        return None
    values = [entry["value"] for entry in entries]
    scale = statistics.fmean([entry["level"] for entry in entries])
    judgement = st.stationarity(
        values, scale=scale, tol_slope_sigma=slope_tolerance
    )
    judgement["slope_tolerance"] = slope_tolerance
    return judgement


def slope_tolerance_for(comparisons, sample_count=None):
    """the slope tolerance, scaled for the screen's many comparisons.

    a two-sigma test applied to twelve independent screens rejects a
    wholly settled bracket most of the time before anything physical has
    happened: the family-wise false-positive rate is one minus
    (0.95) raised to the twelfth, about half. the screen is therefore
    BONFERRONI-corrected, which spends the family's error budget evenly
    across the comparisons.

    the comparison is against a STUDENT'S T quantile rather than a
    normal one, because a fitted slope follows a t law and not a normal
    one. on a fourteen field window the two differ sharply: the normal
    quantile is 2.87 while the t quantile is 3.53. the normal value
    made the screen easier to trip than the mathematics permits, which
    favoured declaring a bracket settled. the degrees of freedom ride
    on the window's own length, so the threshold tightens as the
    window narrows and relaxes as it widens.

    what this widening DOES NOT buy is a guarantee. it pays for the
    twelve screens only. the residuals are serially dependent beyond
    the lag-one term, so the family's real error rate is wider
    than the arithmetic here implies. the rate is measured on surrogate
    records instead of assumed, and the measurement is the number the
    docs quote.
    """
    if comparisons <= 1:
        return sig.bonferroni_threshold(
            1, FAMILY_ALPHA, sample_count or settle_mod.MAX_WINDOW
        )
    return sig.bonferroni_threshold(
        comparisons, FAMILY_ALPHA, sample_count
    )


def assess(rows, core_keys=CORE_KEYS, secondary_keys=SECONDARY_KEYS,
           reference=REFERENCE_TAG):
    """the bracket's settlement, on what each leg adds to the reference."""
    tags = [row["tag"] for row in rows]
    if reference not in tags:
        return {"qualified": False, "reason": "the adiabatic reference is missing"}

    cooled = [tag for tag in tags if tag != reference]
    if not cooled:
        return {"qualified": False, "reason": "no leg has a counterpart to compare"}

    all_keys = tuple(core_keys) + tuple(secondary_keys)
    by_field = series_by_field(rows, all_keys)
    dep = departures_from(by_field, reference, cooled, all_keys)

    for tag in cooled:
        for key in all_keys:
            if key not in dep.get(tag, {}):
                return {
                    "qualified": False,
                    "reason": "%s has no %s series to compare against the reference"
                               % (tag, key),
                }

    shared = None
    for tag in cooled:
        for key in all_keys:
            indices = {entry["index"] for entry in dep[tag][key]}
            shared = indices if shared is None else shared & indices
    if not shared:
        return {"qualified": False, "reason": "the legs share no comparable field"}

    ordered = {
        (tag, key): [e for e in dep[tag][key] if e["index"] in shared]
        for tag in cooled
        for key in all_keys
    }
    unordered = sorted(shared)
    shortest = min(len(series) for series in ordered.values())

    # the screen examines every cooled leg against the reference on every
    # core observable, so the tolerance carries that many comparisons.
    comparisons = max(1, len(cooled) * len(core_keys))
    tolerance = slope_tolerance_for(comparisons)

    # the verdict is made on a DECLARED final horizon. choosing the
    # widest qualifying suffix after seeing the levels lets the data
    # pick its own witness, and the threshold paid for only the twelve
    # declared screens, not that search. the fixed horizon makes the
    # family exactly the screens below and leaves no data-dependent
    # choice to inflate the certificate.
    if shortest < SETTLEMENT_WINDOW:
        return {
            "qualified": False,
            "reason": (
                "the shared record has %d fields, fewer than the declared "
                "%d-field settlement horizon"
                % (shortest, SETTLEMENT_WINDOW)
            ),
        }
    start = shortest - SETTLEMENT_WINDOW
    end = shortest

    # the readings over the chosen window carry that window's own
    # threshold, so a reader can see the number the screen was held to.
    tolerance = slope_tolerance_for(comparisons, end - start)

    core_readings = {}
    secondary_readings = {}
    core_failures = []
    for tag in cooled:
        for key in core_keys:
            judgement = judge(ordered[(tag, key)][start:end], tolerance)
            if judgement is None:
                continue
            core_readings[(tag, key)] = judgement
            if not judgement["qualifies"]:
                core_failures.append("%s/%s: %s" % (tag, key, judgement["reason"]))
        for key in secondary_keys:
            judgement = judge(ordered[(tag, key)][start:end], tolerance)
            if judgement is not None:
                secondary_readings[(tag, key)] = judgement

    # each leg's own LEVEL over the same window. the departure says what
    # the wall changed; the level says where the leg stood, and the
    # gate's choked relation needs the level, not the difference.
    levels = {}
    for tag in tags:
        per_key = {}
        for key in all_keys:
            values = []
            series = (by_field.get(key) or {}).get(tag) or {}
            for index in unordered[start:end]:
                if index in series:
                    values.append(series[index])
            if not values:
                continue
            mean = statistics.fmean(values)
            if len(values) > 1:
                stderr = statistics.stdev(values) / math.sqrt(len(values))
            else:
                stderr = 0.0
            per_key[key] = {
                "mean": mean,
                "standard_error": stderr,
                "effective_n": len(values),
            }
        levels[tag] = per_key

    # the fields the judged window actually covers, in march order.
    window_fields = []
    for tag in cooled:
        for entry in ordered[(tag, core_keys[0])][start:end]:
            window_fields.append(entry["index"])
    window_fields = sorted(set(window_fields))

    out = {
        "qualified": not core_failures,
        "reference": reference,
        "window_start": start,
        "window_end": end,
        "window_length": end - start,
        "window_fields": window_fields,
        "core_keys": list(core_keys),
        "secondary_keys": list(secondary_keys),
        "core": {},
        "secondary": {},
        "levels": levels,
        "slope_tolerance": tolerance,
        "comparisons": comparisons,
        "reason": "; ".join(core_failures),
    }
    for (tag, key), judgement in sorted(core_readings.items()):
        out["core"].setdefault(tag, {})[key] = judgement
    for (tag, key), judgement in sorted(secondary_readings.items()):
        out["secondary"].setdefault(tag, {})[key] = judgement

    unsettled = []
    for tag in cooled:
        for key in secondary_keys:
            judgement = secondary_readings.get((tag, key))
            if judgement is not None and not judgement["qualifies"]:
                unsettled.append("%s/%s: %s" % (tag, key, judgement["reason"]))
    out["secondary_unsettled"] = unsettled
    return out

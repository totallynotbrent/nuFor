"""the settled-window measurement for the bracket's bulk observables.

WHY A SEPARATE STEP

the gate quotes a mean for each observable: the choked mass flow, the
chamber's total pressure and temperature, and the near-wall gas
temperature. that mean is only worth quoting if the series it came
from had settled. this module performs the settlement judgement and
hands back the window it chose, so the gate can say WHICH part of the
run each number describes.

WHAT IT DOES

for each observable it collects the values across the saved fields,
measures the near-wall gas temperature alongside, and asks the
settling criterion whether the tail of the series qualifies as a
stationary window. a qualifier means: the trend across the window is
small against the noise around it, the trend biases the window's mean
by less than the noise, and the trend is statistically
indistinguishable from flat.

once a qualifying window is found, the observable's mean and its
standard error are computed over THAT window, with the effective
sample count correcting for the correlation between successive saved
fields.

WHAT COUNTS AS AN OBSERVABLE HERE

each observable arrives as a list of floats, oldest first. the
stationarity criterion is dimensionless in the observable's own units,
so it applies to a pressure, a mass flow and a temperature alike, with
no per-observable tuning.
"""

import math
import statistics

import stationarity as st

# the window's length, in saved fields.
#
# the bracket saves a field every 250 steps. the flow's transient lasts
# roughly the first 1200 steps, which is about five fields, so a window
# of six fields is the shortest that clears the transient with room to
# spare and still carries ten points. ten points give a trend estimate
# worth trusting and a variance estimate worth quoting, and the
# effective-sample correction handles the rest.
DEFAULT_WINDOW = 10

# the smallest number of fields a window may carry.
#
# a shorter window cannot support a drift-versus-noise comparison, and
# the criterion's own minimum agrees at six. ten is used here because
# the bracket's settled portion needs about five fields to clear the
# transient, so a six-field window would leave only one settled field
# and could not be distinguished from the transient itself.
MIN_WINDOW = 10


def collect(observations, key):
    """the values of one observable across the observations, oldest first."""
    values = []
    for row in observations:
        value = row.get(key)
        if value is None:
            continue
        values.append(float(value))
    return values


def settle(values, window=DEFAULT_WINDOW, **criteria):
    """the settled window of a series, and what it certifies.

    returns a dictionary with the window's indices, its mean, its
    standard error, and the criterion's readings. when no window
    qualifies, the window is empty and `qualified` is false, with the
    reason recorded.
    """
    out = {
        "qualified": False,
        "window_start": None,
        "window_end": None,
        "window_length": 0,
        "mean": None,
        "standard_error": None,
        "effective_n": None,
        "noise_fraction": None,
        "drift_fraction": None,
        "slope_sigma": None,
        "reason": "",
    }

    n = len(values)
    if n < MIN_WINDOW:
        out["reason"] = (
            f"only {n} saved fields, below the {MIN_WINDOW} a settled "
            "window needs"
        )
        return out

    # the candidate windows are contiguous tails of the series. the
    # tail is where the settled state lives, and the window is grown
    # backwards from the last field so the newest data always counts.
    candidates = []
    for length in range(MIN_WINDOW, n + 1):
        candidates.append(length)

    best = None
    for length in reversed(candidates):
        start = n - length
        window_values = values[start:n]
        judgement = st.stationarity(window_values, **criteria)
        if judgement["qualifies"]:
            best = (start, n, judgement)
            break

    if best is None:
        # no tail qualified. judge the fullest permitted window so the
        # reason names the worst offender rather than being blank.
        start = n - MIN_WINDOW
        judgement = st.stationarity(values[start:n], **criteria)
        out["reason"] = judgement["reason"] or "no window settled"
        out["window_start"] = start
        out["window_end"] = n
        out["window_length"] = n - start
        out["mean"] = judgement["mean"]
        out["standard_error"] = judgement["standard_error"]
        out["noise_fraction"] = judgement["noise_fraction"]
        out["drift_fraction"] = judgement["drift_fraction"]
        out["slope_sigma"] = judgement["slope_sigma"]
        out["effective_n"] = judgement["effective_n"]
        return out

    start, end, judgement = best
    out.update(
        {
            "qualified": True,
            "window_start": start,
            "window_end": end,
            "window_length": end - start,
            "mean": judgement["mean"],
            "standard_error": judgement["standard_error"],
            "effective_n": judgement["effective_n"],
            "noise_fraction": judgement["noise_fraction"],
            "drift_fraction": judgement["drift_fraction"],
            "slope_sigma": judgement["slope_sigma"],
        }
    )
    return out


def measure_all(observations, keys, **criteria):
    """the settled window of every observable in a set.

    the observables are settled JOINTLY, not one at a time: the window
    that qualifies for all of them is the one the gate quotes, so the
    numbers describe the same stretch of the run. choosing each
    observable's window separately would let the gate compare a mass
    flow from one stretch with a pressure from another, which is not a
    comparison.
    """
    series = {key: collect(observations, key) for key in keys}
    usable = {key: values for key, values in series.items() if values}
    if not usable:
        return {}

    shortest = min(len(values) for values in usable.values())

    # the joint window: grow it backwards from the last field until
    # EVERY observable qualifies over it.
    chosen = None
    for length in range(shortest, MIN_WINDOW - 1, -1):
        start = shortest - length
        judgements = {}
        ok = True
        for key, values in usable.items():
            judgement = st.stationarity(values[start:], **criteria)
            judgements[key] = judgement
            if not judgement["qualifies"]:
                ok = False
                break
        if ok:
            chosen = (start, length, judgements)
            break

    if chosen is None:
        # nothing satisfied every observable at once. fall back on the
        # most permissive reading available and report which observable
        # objected, so the gate can say the window is not certified.
        start = max(0, shortest - MAX_WINDOW)
        judgements = {
            key: st.stationarity(values[start:], **criteria)
            for key, values in usable.items()
        }
        blockers = [
            f"{key}: {judgement['reason']}"
            for key, judgement in judgements.items()
            if not judgement["qualifies"]
        ]
        out = {
            "qualified": False,
            "window_start": start,
            "window_end": shortest,
            "window_length": shortest - start,
            "observables": {},
            "reason": "; ".join(blockers),
        }
        for key, judgement in judgements.items():
            out["observables"][key] = {
                "mean": judgement["mean"],
                "standard_error": judgement["standard_error"],
                "effective_n": judgement["effective_n"],
                "noise_fraction": judgement["noise_fraction"],
                "drift_fraction": judgement["drift_fraction"],
                "slope_sigma": judgement["slope_sigma"],
                "qualifies": judgement["qualifies"],
            }
        return out

    start, length, judgements = chosen
    out = {
        "qualified": True,
        "window_start": start,
        "window_end": shortest,
        "window_length": length,
        "observables": {},
        "reason": "",
    }
    for key, judgement in judgements.items():
        out["observables"][key] = {
            "mean": judgement["mean"],
            "standard_error": judgement["standard_error"],
            "effective_n": judgement["effective_n"],
            "noise_fraction": judgement["noise_fraction"],
            "drift_fraction": judgement["drift_fraction"],
            "slope_sigma": judgement["slope_sigma"],
            "qualifies": judgement["qualifies"],
        }
    return out


# the widest window considered when no shorter one qualifies.
MAX_WINDOW = 30


def paired_contrast(reference, other, key, sem_lookup=None):
    """the difference between two legs' observables, with its significance.

    the two legs are marched from the same initial condition on the
    same grid, and their fields are saved at the same steps, so the
    pairing is by FIELD INDEX and the difference between the legs at
    each field cancels everything the wall did not change.

    the contrast's significance uses the difference series' OWN spread,
    not the two legs' spreads added in quadrature. the legs share their
    transients almost entirely, so the difference is far quieter than
    either leg alone, and treating them as independent throws away
    most of the sensitivity. the effective sample count corrects for
    the residual correlation between successive differences.
    """
    ref_values = reference.get(key + "_series") or []
    oth_values = other.get(key + "_series") or []
    n = min(len(ref_values), len(oth_values))
    if n < 3:
        return None

    differences = [oth_values[i] - ref_values[i] for i in range(n)]

    # restrict to the window the legs were quoted over, so the contrast
    # describes the same stretch of the run as the levels themselves.
    start = other.get("_window_start", 0)
    if start < n:
        differences = differences[start:n]
        n = len(differences)
    if n < 3:
        return None

    mean_difference = statistics.fmean(differences)
    deviations = [d - mean_difference for d in differences]
    variance = sum(d * d for d in deviations) / (n - 1)
    sigma = math.sqrt(max(variance, 0.0))
    effective_n = st.bartlett_neffective(deviations)
    sem = sigma / math.sqrt(effective_n) if effective_n > 0 else float("inf")

    reference_mean = statistics.fmean(ref_values[start:])
    scale = abs(reference_mean) if abs(reference_mean) > 0.0 else 1.0

    return {
        "delta": mean_difference,
        "delta_pct": 100.0 * mean_difference / scale,
        "sd": sigma,
        "effective_n": effective_n,
        "standard_error": sem,
        "sigma": mean_difference / sem if sem > 0.0 else float("inf"),
        "reference_mean": reference_mean,
        "n": n,
    }

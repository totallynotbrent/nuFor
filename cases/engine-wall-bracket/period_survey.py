"""the legacy period detector, kept only to REPORT the absence of a period.

WHY THIS STILL EXISTS

an earlier version of this gate inferred the flow's period from the
saved fields and averaged over whole periods. measurement refuted the
underlying premise: the fields saved at steps 250 and 3750 have
different checksums, and the correlation at every tested lag is about
zero. the inviscid predecessor in `cases/chamber-throat` shows the same
damped rise-then-decay with decaying oscillatory residuals, so the
viscosity did not suppress a real oscillation; the premise was
borrowed from the wrong case.

the settled-window criterion in `stationarity.py` replaces the
averaging. this module remains because the ABSENCE of a period is
itself worth stating in the gate's report: a reader should see that
the bracket was checked for periodicity and that the flow was found
NOT to be periodic, rather than seeing the question silently dropped.

nothing here feeds a contract. `detect_period_across_legs` and its
helpers produce the report's descriptive line only.
"""

import glob
import hashlib
import math
import os

# the largest lag the detector will consider.
#
# a field is saved every 250 steps and a periodic orbit in this flow
# would have spanned tens of fields, so a search out to ninety-six
# fields covers a period of twenty-four thousand steps: four times
# the transient's own length.
MAX_LAG = 96

# the smallest overlap a lag may keep and still be believed.
#
# a correlation computed from a handful of overlapping pairs is a
# coincidence dressed as a signal: the tails of the search reach
# lags whose overlap is a few fields wide, and there the numbers swing
# hard for no reason in the flow. a lag has to keep at least this many
# pairs, and at least this share of the series, before its correlation
# is allowed to speak.
MIN_OVERLAP = 8
MIN_OVERLAP_SHARE = 0.25


def digest(path):
    """the field's checksum."""
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            hasher.update(block)
    return hasher.hexdigest()


def signal_correlation(values, lag):
    """the correlation between a series and itself, displaced by a lag.

    the series is detrended by subtracting its own linear fit before
    the correlation is taken, so a series that merely drifts does not
    read as an oscillator.
    """
    n = len(values)
    if lag <= 0 or n - lag < 3:
        return 0.0
    mean_x = sum(range(n)) / n
    mean_y = sum(values) / n
    sxx = sum((i - mean_x) ** 2 for i in range(n))
    syy = sum((v - mean_y) ** 2 for v in values)
    if sxx <= 0.0 or syy <= 0.0:
        return 0.0
    sxy = sum((i - mean_x) * (values[i] - mean_y) for i in range(n))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    residual = [values[i] - (intercept + slope * i) for i in range(n)]

    a = residual[:-lag]
    b = residual[lag:]
    ma = sum(a) / len(a)
    mb = sum(b) / len(b)
    va = sum((x - ma) ** 2 for x in a)
    vb = sum((x - mb) ** 2 for x in b)
    if va <= 0.0 or vb <= 0.0:
        return 0.0
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    return num / math.sqrt(va * vb)



def detect_period(digests):
    """the flow's period in saved fields, from the fields' checksums.

    a strict PERIODIC orbit repeats a field exactly, so two saved
    fields with the same checksum mark one whole period apart. a
    settled flow repeats nothing, and this returns None.
    """
    seen = {}
    for index, value in enumerate(digests):
        if value in seen:
            gap = index - seen[value]
            if gap > 0:
                return gap, seen[value], index
        seen[value] = index
    return None, None, None


def detect_period_by_signal(values):
    """the period suggested by an observable's own autocorrelation.

    this is the weaker of the two tests, and it is prone to suggesting
    spurious lags from a short series whose trend mimics a cycle. it
    is kept for the report, where the value of a LOW best correlation
    is the point.

    only lags that keep enough overlapping pairs are considered. the
    search reaches far past the series' own length, and there the
    surviving pairs are numbered in single digits, which makes the
    correlation swing hard for reasons that have nothing to do with the
    flow.
    """
    if len(values) < 4:
        return None, 0.0
    count = len(values)
    floor = max(MIN_OVERLAP, int(count * MIN_OVERLAP_SHARE))
    best_lag = None
    best_score = 0.0
    for lag in range(1, min(MAX_LAG, count - 3)):
        if count - lag < floor:
            break
        score = signal_correlation(values, lag)
        if score > best_score:
            best_lag = lag
            best_score = score
    if best_lag is None or best_score <= 0.0:
        return None, best_score
    return best_lag, best_score


def survey(per_leg):
    """whether the bracket's legs are periodic, for the report only.

    the legs march the same chamber, so a genuine period would be
    shared. the survey reports, for each leg, the checksum repeat and
    the strongest autocorrelation, so the gate's report can say that
    the flows were checked and found not to oscillate.

    `per_leg` maps a leg's name to its saved fields' checksums and the
    mass flow read from them. the caller supplies both, so the survey
    reads the same way over the marched field and over the record.
    """
    out = {}
    for tag, fields in per_leg.items():
        digests = fields.get("digests") or []
        values = fields.get("values") or []
        if not digests and not values:
            out[tag] = {"checked": False, "why": "no saved fields"}
            continue
        gap, _, _ = detect_period(digests)
        best_lag, best_score = detect_period_by_signal(values)
        out[tag] = {
            "checked": True,
            "fields": len(digests) or len(values),
            "checksum_period": gap,
            "best_autocorrelation_lag": best_lag,
            "best_autocorrelation": best_score,
            "periodic": gap is not None,
        }
    return out

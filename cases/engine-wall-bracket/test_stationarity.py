"""check the settling criterion against series whose behaviour is known.

these are behavioural checks, written against what the criterion must
judge, not snapshots of its current numbers. each case supplies a series
whose character is set by construction: a settled wander, a steady
trend, a sustained oscillation, a single outlier, or a stepped change.
the criterion's answer on each is determined by that construction, so a
future change to its internals cannot quietly redefine it.

run from this directory:

    python3.12 test_stationarity.py
"""

import math
import random
import sys

import stationarity as st


FAILURES = []


def maths_is_finite(value):
    return isinstance(value, float) and math.isfinite(value)


def check(condition, message):
    if not condition:
        FAILURES.append(message)


def expect_qualifies(name, series, **criteria):
    judgement = st.stationarity(series, **criteria)
    check(
        judgement["qualifies"],
        f"{name}: expected the window to qualify, but it did not ({judgement['reason']})",
    )
    return judgement


def expect_rejects(name, series, **criteria):
    judgement = st.stationarity(series, **criteria)
    check(
        not judgement["qualifies"],
        f"{name}: expected the window to be rejected, but it qualified",
    )
    return judgement


def main():
    rng = random.Random(20261101)

    # a settled wander: noise about a level, with no drift. the
    # criterion must accept it.
    settled = [100.0 + rng.gauss(0.0, 0.5) for _ in range(40)]
    judgement = expect_qualifies("a settled wander", settled)
    check(
        judgement["noise_fraction"] > 0.0,
        "a settled wander should report a noise fraction",
    )

    # the same wander with a mild drift: two thousandths of the level
    # per sample. over forty samples that is eight percent of the mean,
    # which the mean cannot hide.
    drifting = [
        100.0 + 0.02 * i + rng.gauss(0.0, 0.5) for i in range(40)
    ]
    expect_rejects("a steadily drifting series", drifting)

    # a sustained oscillation: the flow's actual character. the
    # criterion cares about a TREND, not about variation, so a bounded
    # wobble is not grounds for rejecting a window. the drift through
    # an integer number of cycles is nearly nil.
    oscillating = [100.0 + 2.0 * math.sin(i * 0.7) for i in range(40)]
    expect_qualifies("a sustained oscillation", oscillating)

    # a ramp with no noise at all: unmistakable, and must be rejected.
    ramp = [float(i) for i in range(40)]
    judgement = expect_rejects("a noiseless ramp", ramp)
    check(
        judgement["drift_fraction"] > 0.0,
        "a ramp's drift fraction should be positive",
    )

    # a single gross outlier among settled points. the robust reading
    # is that the trend is still nil, so the window's MEAN may be
    # skewed, but the criterion judges settling rather than outliers.
    # the important property is that the outlier does not manufacture a
    # significant TREND.
    spiked = [100.0 + rng.gauss(0.0, 0.3) for _ in range(39)] + [180.0]
    judgement = st.stationarity(spiked)
    check(
        judgement["drift_fraction"] < 0.02,
        "one outlier should not read as a drift",
    )

    # a step halfway through: a settled level, a jump, then another
    # settled level. the window straddles the discontinuity, so a
    # linear fit across it reads a spurious trend and the window must
    # be rejected rather than averaged across the step.
    stepped = (
        [100.0 + rng.gauss(0.0, 0.3) for _ in range(20)]
        + [120.0 + rng.gauss(0.0, 0.3) for _ in range(20)]
    )
    expect_rejects("a window straddling a step", stepped)

    # a perfectly flat series has no noise and no drift. it must
    # qualify, and the divisions by zero must be handled rather than
    # raising.
    flat = [50.0] * 30
    judgement = expect_qualifies("a perfectly flat series", flat)
    check(
        maths_is_finite(judgement["standard_error"]),
        "a flat series should report a finite standard error",
    )

    # the window search takes the LONGEST qualifying window, because a
    # longer window carries less random error on its mean.
    mixed = [100.0 + rng.gauss(0.0, 0.5) for _ in range(30)] + [
        100.0 + 0.3 * i + rng.gauss(0.0, 0.5) for i in range(20)
    ]
    length, judgement = st.best_window(mixed)
    check(
        judgement["qualifies"],
        "the search should find the settled stretch at the start of the series",
    )
    check(
        0 < length <= 30,
        f"the chosen window should lie inside the settled stretch, got {length}",
    )
    check(
        judgement["length"] == length,
        "the returned judgement should carry the chosen length",
    )

    # a series with no settled stretch at all: the search reports the
    # failure instead of silently returning a bad window.
    hopeless = [100.0 + 0.5 * i + rng.gauss(0.0, 0.5) for i in range(40)]
    length, judgement = st.best_window(hopeless)
    check(length == 0, "a wholly drifting series should yield no window")
    check(
        bool(judgement["reason"]),
        "a failed search should carry a reason",
    )

    # a series shorter than the minimum cannot be judged, and must say
    # so rather than judging it anyway.
    judgement = st.stationarity([1.0, 2.0, 3.0])
    check(not judgement["qualifies"], "a tiny series cannot qualify")
    check("below" in judgement["reason"], "a tiny series should say why")

    # the effective sample count must fall when the samples correlate
    # and stay near the raw count when they do not.
    independent = [rng.gauss(0.0, 1.0) for _ in range(60)]
    eff_independent = st.bartlett_neffective(independent)
    check(
        eff_independent > 30.0,
        f"independent samples should keep most of their count, got {eff_independent:.1f}",
    )
    correlated = []
    value = 0.0
    for _ in range(60):
        value = 0.95 * value + rng.gauss(0.0, 0.3)
        correlated.append(value)
    eff_correlated = st.bartlett_neffective(correlated)
    check(
        eff_correlated < eff_independent,
        f"correlated samples should lose count, got {eff_correlated:.1f} vs {eff_independent:.1f}",
    )

    if FAILURES:
        print(f"FAILED: {len(FAILURES)} check(s)")
        for failure in FAILURES:
            print(f"  - {failure}")
        return 1
    print("stationarity: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

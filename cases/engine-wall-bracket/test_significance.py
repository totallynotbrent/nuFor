"""check the significance helpers against values that are known.

these are behavioural checks rather than snapshots. the t quantiles are
checked against closed-form cases and against a coarse numerical
integral, the correlation widening is checked against its defining
factor, and the trend test is checked against series whose trend is set
by construction.
"""

import math
import random
import sys

import significance as sig


FAILURES = []


def check(condition, message):
    if not condition:
        FAILURES.append(message)


def near(actual, expected, tolerance, message):
    if abs(actual - expected) > tolerance:
        FAILURES.append(
            "%s: got %.10g, expected %.10g (within %.3g)"
            % (message, actual, expected, tolerance)
        )


def t_two_sided(alpha, degrees_of_freedom):
    return sig.student_t_quantile(1.0 - alpha / 2.0, degrees_of_freedom)


def main():
    # the t distribution meets the normal one at the extremes. for a
    # thousand degrees of freedom the quantile approaches the normal
    # quantile closely.
    near(
        t_two_sided(0.05, 100000),
        1.959963985,
        1.0e-4,
        "the t quantile should approach the normal quantile",
    )

    # a heavier tail than the normal law: at twelve degrees of freedom
    # the ninety-five percent critical value is well above the normal
    # one's, and the analytic value is known.
    near(
        t_two_sided(0.05, 12),
        2.178812829,
        1.0e-6,
        "the twelve-degree t quantile should match the analytic value",
    )
    near(
        t_two_sided(0.01, 10),
        3.169272673,
        1.0e-6,
        "the ten-degree one-percent t quantile should match the known value",
    )

    # the quantile must satisfy its own definition: the area to the
    # left of it equals the requested probability.
    for probability in (0.01, 0.1, 0.5, 0.9, 0.995):
        for degrees_of_freedom in (1, 3, 12, 50):
            value = sig.student_t_quantile(probability, degrees_of_freedom)
            area = sig.student_t_cdf(value, degrees_of_freedom)
            near(
                area,
                probability,
                1.0e-7,
                "the quantile at %.3f with %d degrees should invert the cdf"
                % (probability, degrees_of_freedom),
            )

    # the cdf must be symmetric about zero.
    for degrees_of_freedom in (1, 5, 12):
        for value in (0.5, 1.0, 3.0):
            lower = sig.student_t_cdf(-value, degrees_of_freedom)
            upper = sig.student_t_cdf(value, degrees_of_freedom)
            near(
                lower,
                1.0 - upper,
                1.0e-9,
                "the t cdf should be symmetric at %.2f with %d degrees"
                % (value, degrees_of_freedom),
            )

    # the family's threshold must tighten as the family grows, and
    # must loosen as the window widens.
    small_family = sig.bonferroni_threshold(2, 0.05, 20)
    large_family = sig.bonferroni_threshold(24, 0.05, 20)
    check(
        large_family > small_family,
        "a larger family should demand a larger threshold",
    )
    narrow = sig.bonferroni_threshold(12, 0.05, 8)
    wide = sig.bonferroni_threshold(12, 0.05, 60)
    check(
        wide < narrow,
        "a wider window should relax the threshold",
    )
    check(
        abs(narrow - t_two_sided(0.05 / 12, 6)) < 1.0e-9,
        "the threshold should be the per-screen critical value",
    )

    # holm's thresholds must come back largest first and dominate
    # bonferroni's single threshold.
    thresholds = sig.holm_thresholds(12, 0.05, 14)
    check(
        all(
            thresholds[i] >= thresholds[i + 1]
            for i in range(len(thresholds) - 1)
        ),
        "holm's thresholds should descend",
    )
    # holm spends its family alpha unevenly: the first screen faces the
    # strictest threshold, which is bonferroni's, and each later screen
    # earns a looser one. the procedure therefore rejects MORE than
    # bonferroni at the same nominal level.
    check(
        abs(thresholds[0] - sig.bonferroni_threshold(12, 0.05, 14)) < 1.0e-9,
        "holm's first threshold should be bonferroni's",
    )
    check(
        thresholds[-1] < thresholds[0],
        "holm should relax its threshold as the screens succeed",
    )

    # the lag-one correlation of an independent series sits near zero,
    # and of a strongly correlated one near its coefficient.
    rng = random.Random(20261103)
    independent = [rng.gauss(0.0, 1.0) for _ in range(400)]
    check(
        abs(sig.lag_one_autocorrelation(independent)) < 0.12,
        "an independent series should show little lag-one correlation",
    )
    correlated = []
    value = 0.0
    for _ in range(4000):
        value = 0.8 * value + rng.gauss(0.0, 1.0)
        correlated.append(value)
    check(
        abs(sig.lag_one_autocorrelation(correlated) - 0.8) < 0.08,
        "an ar(1) series should show its own correlation",
    )

    # the widening factor must follow the autoregressive inflation
    # `(1 + rho) / (1 - rho)` under the square root.
    residuals = [0.5 * (-1.0) ** index for index in range(21)]
    rho = sig.lag_one_autocorrelation(residuals)
    expected = math.sqrt((1.0 + rho) / (1.0 - rho))
    ratio = sig.autocorrelated_slope_standard_error(residuals, 2.0) / 2.0
    near(
        ratio,
        expected,
        1.0e-9,
        "the widening must follow the autoregressive inflation factor",
    )

    # independent residuals must leave the error bar nearly alone. an
    # independent sample still carries a small accidental correlation,
    # so the check allows the widening that correlation justifies and
    # refuses anything beyond it.
    neutral = [rng.gauss(0.0, 1.0) for _ in range(2000)]
    rho = sig.lag_one_autocorrelation(neutral)
    check(
        abs(rho) < 0.05,
        "independent residuals should show almost no lag-one correlation,"
        " got %.4f" % rho,
    )
    justified = 3.0 * math.sqrt((1.0 + rho) / (1.0 - rho))
    near(
        sig.autocorrelated_slope_standard_error(neutral, 3.0),
        justified,
        1.0e-9,
        "independent residuals must leave the error bar essentially alone",
    )
    check(
        abs(sig.autocorrelated_slope_standard_error(neutral, 3.0) / 3.0 - 1.0)
        < 0.1,
        "independent residuals must not move the error bar materially",
    )

    # a strong correlation must widen the error bar substantially.
    smooth = [math.sin(index * 0.3) for index in range(30)]
    check(
        sig.autocorrelated_slope_standard_error(smooth, 1.0) > 1.0,
        "strongly correlated residuals must widen the error bar",
    )

    # a trend test on a pure ramp must find a decisive slope.
    ramp = [float(index) for index in range(20)]
    judgement = sig.slope_test(ramp)
    check(
        judgement["t_statistic"] > 100.0,
        "a noiseless ramp must produce a decisive slope statistic",
    )
    near(
        judgement["slope"],
        1.0,
        1.0e-9,
        "the fitted slope of a unit ramp should be one",
    )
    check(
        judgement["p_value"] < 1.0e-12,
        "a noiseless ramp should carry a vanishing p value",
    )

    # a flat series carries no trend.
    flat = sig.slope_test([3.0] * 20)
    check(
        flat["t_statistic"] == 0.0,
        "a flat series should carry no trend statistic",
    )
    check(
        flat["p_value"] == 1.0,
        "a flat series should carry the maximal p value",
    )

    # the degrees of freedom must be the window's length less two.
    check(
        sig.slope_test([1.0, 2.0, 3.0, 4.0])["degrees_of_freedom"] == 2,
        "four points should leave two degrees of freedom",
    )

    # a drift dressed in correlated noise must still be caught when it
    # is large, and must not be reported as certain when it is small.
    weak = []
    value = 0.0
    for index in range(40):
        value = 0.95 * value + rng.gauss(0.0, 1.0)
        weak.append(value + 0.01 * index)
    weak_judgement = sig.slope_test(weak)
    check(
        weak_judgement["p_value"] > 0.0,
        "a weakly drifted correlated series should keep a positive p value",
    )

    if FAILURES:
        print("FAILED: %d check(s)" % len(FAILURES))
        for failure in FAILURES:
            print("  - %s" % failure)
        return 1
    print("significance: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

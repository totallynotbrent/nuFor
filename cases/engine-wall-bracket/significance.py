"""significance testing for the bracket's settlement screens.

why this module exists

the bracket judges each wall's departure from the adiabatic reference
with a slope test. two faults sat in that test. it compared the slope
against a normal quantile, where a regression slope follows a student's
t distribution, and it used a slope standard error that ignored the
correlation between neighbouring saved fields. both faults favoured the
bracket: a smaller threshold than the data warrants, and an error bar
too narrow because correlated samples were counted as independent.

what this module provides

- a student's t quantile computed from the standard library alone;
- a slope standard error widened for the residuals' lag-one correlation,
  following the bartlett correction the bracket already uses elsewhere;
- the family-wise threshold for a screen of many comparisons.

the functions here are stateless and take plain lists, so the bracket's
scripts gain the correction without gaining a dependency.
"""

import math

# the smallest degree of freedom the t distribution accepts.
_MIN_DEGREES_OF_FREEDOM = 1

# the width of the quantile search interval.
_QUANTILE_EPSILON = 1.0e-10

# the extreme values the quantile search brackets.
_QUANTILE_LOW = -1.0e3
_QUANTILE_HIGH = 1.0e3

# the iteration cap for the continued fraction expansion.
_FRACTION_STEPS = 300

# the floor that keeps the continued fraction division safe.
_TINY = 1.0e-300


def _beta_continued_fraction(a, b, x):
    """the continued fraction expansion of the incomplete beta function."""
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < _TINY:
        d = _TINY
    d = 1.0 / d
    h = d
    for m in range(1, _FRACTION_STEPS + 1):
        doubled = 2 * m
        term = m * (b - m) * x / ((qam + doubled) * (a + doubled))
        d = 1.0 + term * d
        if abs(d) < _TINY:
            d = _TINY
        c = 1.0 + term / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        h *= d * c
        term = -(a + m) * (qab + m) * x / ((a + doubled) * (qap + doubled))
        d = 1.0 + term * d
        if abs(d) < _TINY:
            d = _TINY
        c = 1.0 + term / c
        if abs(c) < _TINY:
            c = _TINY
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 3.0e-16:
            break
    return h


def _regularised_beta(a, b, x):
    """the regularised incomplete beta function."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_front = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log1p(-x)
    )
    front = math.exp(log_front)
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_continued_fraction(a, b, x) / a
    return 1.0 - front * _beta_continued_fraction(b, a, 1.0 - x) / b


def student_t_cdf(statistic, degrees_of_freedom):
    """the area to the left of a value under the student's t law."""
    if degrees_of_freedom < _MIN_DEGREES_OF_FREEDOM:
        raise ValueError("degrees of freedom must be at least one")
    if statistic == 0.0:
        return 0.5
    ratio = degrees_of_freedom / (
        degrees_of_freedom + statistic * statistic
    )
    tail = 0.5 * _regularised_beta(
        degrees_of_freedom / 2.0, 0.5, ratio
    )
    if statistic > 0.0:
        return 1.0 - tail
    return tail


def student_t_quantile(probability, degrees_of_freedom):
    """the value with `probability` of the student's t law below it.

    the two-sided critical value for a level `alpha` is
    `student_t_quantile(1.0 - alpha / 2.0, degrees_of_freedom)`.
    """
    if not 0.0 < probability < 1.0:
        raise ValueError("the probability must lie between zero and one")
    if degrees_of_freedom < _MIN_DEGREES_OF_FREEDOM:
        raise ValueError("degrees of freedom must be at least one")
    low = _QUANTILE_LOW
    high = _QUANTILE_HIGH
    for _ in range(200):
        middle = 0.5 * (low + high)
        if student_t_cdf(middle, degrees_of_freedom) < probability:
            low = middle
        else:
            high = middle
        if high - low < _QUANTILE_EPSILON:
            break
    return 0.5 * (low + high)


def lag_one_autocorrelation(residuals):
    """the lag-one autocorrelation of a residual series."""
    count = len(residuals)
    if count < 3:
        return 0.0
    mean = sum(residuals) / count
    centred = [residual - mean for residual in residuals]
    denominator = sum(value * value for value in centred)
    if denominator <= 0.0:
        return 0.0
    numerator = sum(
        centred[index] * centred[index + 1] for index in range(count - 1)
    )
    return numerator / denominator


def autocorrelated_slope_standard_error(residuals, raw_standard_error):
    """the slope's standard error, widened for lag-one correlation.

    for a first order autoregressive error the variance of an ordinary
    least squares slope is inflated by `(1 + rho) / (1 - rho)`, where
    `rho` is the residuals' lag-one autocorrelation. an uncorrelated
    series leaves the error untouched.
    """
    rho = lag_one_autocorrelation(residuals)
    if rho >= 0.999:
        return float("inf")
    factor = (1.0 + rho) / (1.0 - rho)
    return raw_standard_error * math.sqrt(factor)


def ols_slope(values):
    """the least squares slope of a series against its own index."""
    count = len(values)
    if count < 2:
        raise ValueError("a fit needs at least two points")
    mean_x = (count - 1) / 2.0
    mean_y = sum(values) / count
    sxx = sum((x - mean_x) ** 2 for x in range(count))
    if sxx <= 0.0:
        raise ValueError("the abscissa must vary")
    sxy = sum(
        (x - mean_x) * (y - mean_y) for x, y in enumerate(values)
    )
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    residuals = [
        y - (intercept + slope * x) for x, y in enumerate(values)
    ]
    return slope, intercept, residuals


def slope_test(values):
    """the significance of a window's trend, with the correlations kept.

    the returned dictionary carries the slope, the residuals' lag-one
    autocorrelation, the raw and corrected standard errors, the t
    statistic and the two sided p value.
    """
    count = len(values)
    if count < 3:
        raise ValueError("a trend test needs at least three points")
    slope, _, residuals = ols_slope(values)
    degrees_of_freedom = count - 2
    residual_sum_squares = sum(r * r for r in residuals)
    sigma = math.sqrt(max(residual_sum_squares / degrees_of_freedom, 0.0))
    sxx = sum((x - (count - 1) / 2.0) ** 2 for x in range(count))
    if sigma > 0.0 and sxx > 0.0:
        raw_standard_error = sigma / math.sqrt(sxx)
    elif slope == 0.0:
        raw_standard_error = 0.0
    else:
        raw_standard_error = float("inf")
    standard_error = autocorrelated_slope_standard_error(
        residuals, raw_standard_error
    )
    if standard_error == 0.0:
        statistic = 0.0 if slope == 0.0 else float("inf")
    elif math.isinf(standard_error):
        statistic = float("inf")
    else:
        statistic = slope / standard_error
    p_value = 2.0 * (
        1.0 - student_t_cdf(abs(statistic), degrees_of_freedom)
    )
    return {
        "n": count,
        "slope": slope,
        "residual_acf1": lag_one_autocorrelation(residuals),
        "slope_se_raw": raw_standard_error,
        "slope_se": standard_error,
        "t_statistic": statistic,
        "degrees_of_freedom": degrees_of_freedom,
        "p_value": min(max(p_value, 0.0), 1.0),
    }


def bonferroni_threshold(
    comparisons, family_alpha=0.05, sample_count=None
):
    """the student's t threshold for a family of screens.

    `comparisons` counts the screens in the family and `sample_count`
    gives the window's length, which sets the degrees of freedom.
    """
    if comparisons <= 1:
        return student_t_quantile(1.0 - family_alpha / 2.0, 1000)
    per_screen_alpha = family_alpha / comparisons
    degrees_of_freedom = 1000
    if sample_count is not None:
        degrees_of_freedom = sample_count - 2
    if degrees_of_freedom < _MIN_DEGREES_OF_FREEDOM:
        degrees_of_freedom = _MIN_DEGREES_OF_FREEDOM
    return student_t_quantile(
        1.0 - per_screen_alpha / 2.0, degrees_of_freedom
    )


def holm_thresholds(
    comparisons, family_alpha=0.05, sample_count=None
):
    """the student's t thresholds of holm's step-down procedure.

    the thresholds come back largest first, so the largest p value faces
    the loosest threshold. holm dominates bonferroni and is the better
    choice when the screens are not equally weighted.
    """
    degrees_of_freedom = 1000
    if sample_count is not None:
        degrees_of_freedom = sample_count - 2
    if degrees_of_freedom < _MIN_DEGREES_OF_FREEDOM:
        degrees_of_freedom = _MIN_DEGREES_OF_FREEDOM
    thresholds = []
    for rank in range(comparisons):
        alpha = family_alpha / (comparisons - rank)
        thresholds.append(
            student_t_quantile(1.0 - alpha / 2.0, degrees_of_freedom)
        )
    return thresholds

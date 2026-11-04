"""the bracket's settling criterion: statistical stationarity.

WHY THIS EXISTS

an earlier version of this gate assumed the marched chamber fell into a
strictly periodic orbit, with the field saved at step 250 byte-identical
to the one saved at step 3750. measurement refuted that on both counts.
the fields at those steps have different checksums, and the correlation
of the throat's mass flow at every tested lag is about zero. the
inviscid predecessor in `cases/chamber-throat` shows the same damped
rise-then-decay with decaying oscillatory residuals, so the bracket's
viscosity did not suppress a real oscillation; the premise was borrowed
from the wrong case.

what the flow actually does is damp toward a statistical stationary
state: the observable wanders about a slowly varying level, with a
scatter of some tenths of a percent and a drift of a few hundredths of
a percent per field. averaging over such a window is legitimate, and
that is what this module decides. a window qualifies only when it is
long enough that the drift across it is small compared with the noise
around it, and when the trend is not distinguishable from zero.

WHAT THE CRITERION IS

three quantities, all measured on one observable's series:

- DRIFT: the slope of a least squares fit across the window, expressed
  in the observable's own units per field, and as a percentage of the
  window's mean per field. this is the bias that averaging cannot undo.

- NOISE: the standard deviation of the residuals about that fit,
  divided by the window's mean. this is the scatter the mean sits in.

- SETTLING: the drift's magnitude against the noise, tested as a
  regression slope. a slope significantly different from zero means the
  window has a trend, and the mean over it is biased.

a window QUALIFIES when the drift's contribution to the mean's
uncertainty is negligible beside the noise, and the slope's confidence
interval excludes any meaningful trend. the qualification is written in
terms of the ratio of drift to noise, so it is dimensionless and
carries across observables and across the bracket's legs.

WHY THE WINDOW MATTERS

the drift accumulates. a window twice as long carries twice the drift
but no more than sqrt(2) times the mean's random error, so a long
window is always quieter. the trade is that a long window spends more
of the run on a leg that may still be moving. the criterion balances
those directly: it asks whether the trend is small enough, over the
window actually chosen, that the mean can be quoted with a credible
uncertainty.

TIME AUTOCORRELATION

the samples are not independent. a field every 250 steps inherits the
previous one's state, so the effective sample count is below the raw
count and the naive standard error UNDERSTATES the uncertainty. this
uses Bartlett's independence test on the residual series to estimate the
effective count, which inflates the error bars where the series is
strongly autocorrelated and leaves them alone where it is not. quoting
a mean without it would overstate the precision of the bracket.

this module deals only in the statistics. the gate's own contracts,
the observables it reads, and the physical checks it asserts live in
`gates.py`.
"""

import math

import significance as sig

# the minimum number of samples a qualifying window needs.
#
# one sample cannot carry a variance and two cannot carry a trend, so
# a window shorter than this cannot be judged at all. the number is
# generous rather than strict: the judge is the drift-to-noise ratio,
# which is meaningless on a handful of noisy points.
MIN_WINDOW = 6


def ols_fit(xs, ys):
    """the least squares slope and intercept of ys against xs.

    returns (slope, intercept, residuals, n). the slope's standard
    error needs the residual variance, so the residuals come back too.
    """
    n = len(xs)
    if n < 2:
        raise ValueError("a fit needs at least two points")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    if sxx <= 0.0:
        raise ValueError("a fit needs the abscissa to vary")
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    residuals = [y - (intercept + slope * x) for x, y in zip(xs, ys)]
    return slope, intercept, residuals, n


def bartlett_neffective(residuals):
    """the effective sample count of a correlated residual series.

    Bartlett's independence test compares the lag-one autocorrelation
    against the value chance alone would give. summing those ratios
    over the lags gives the factor by which the true sample count
    exceeds the apparent one; the effective count is the apparent count
    divided by it. a series of independent samples gives a factor near
    one, and a smoothly drifting series gives a large one.
    """
    n = len(residuals)
    if n < 3:
        return float(n)
    var = sum((r - sum(residuals) / n) ** 2 for r in residuals) / (n - 1)
    if var <= 0.0:
        return float(n)
    factor = 0.0
    for lag in range(1, n):
        cov = sum(
            residuals[i] * residuals[i + lag] for i in range(n - lag)
        ) / (n - lag)
        rho = cov / var
        if rho <= 0.0:
            break
        factor += (1.0 - lag / n) * rho
        if factor > n:
            factor = float(n)
            break
    factor = max(factor, 1.0)
    return n / factor


def stationarity(
    values, tol_drift_fraction=0.02, tol_slope_sigma=2.0, scale=None
):
    """judge whether a window of a series has settled.

    the observable is `values`, oldest first. the judgement combines
    three readings of the same window:

    - the drift, which is the fitted slope per sample, and which biases
      the window's mean;
    - the noise, which is the scatter about the fit, and which the mean
      sits inside;
    - the slope's own significance, so a window that trends however
      gently is rejected rather than quietly averaged.

    `tol_drift_fraction` bounds the drift as a fraction of the mean per
    sample: at 0.02, a window whose level moves by two percent of its
    own mean across its whole length is already suspect. the criterion
    actually applied is stricter, because it asks the drift's bias on
    the mean to be small beside the noise.

    `tol_slope_sigma` is the number of standard errors the slope may
    sit from zero before the window counts as trending.

    `scale` is the magnitude the drift and the noise are measured
    against. it defaults to the series' own mean, which suits a level.
    a difference series has a mean near zero, so passing the level it
    was differenced from judges the difference fairly instead of dividing
    by almost nothing.

    returns a dictionary with the readings and a `qualifies` flag. the
    caller decides what to do with a window that fails.
    """
    n = len(values)
    out = {
        "n": n,
        "mean": None,
        "noise_fraction": None,
        "drift_fraction": None,
        "drift_to_noise": None,
        "slope_sigma": None,
        "effective_n": None,
        "standard_error": None,
        "drift_bias_fraction": None,
        "qualifies": False,
        "reason": "",
    }
    if n < MIN_WINDOW:
        out["reason"] = f"only {n} samples, below the {MIN_WINDOW} needed"
        return out

    xs = list(range(n))
    mean = sum(values) / n
    slope, _, residuals, _ = ols_fit(xs, values)

    # the scatter about the fit, in the observable's own units.
    resid_var = sum(r * r for r in residuals) / (n - 2)
    sigma = math.sqrt(max(resid_var, 0.0))

    eff_n = bartlett_neffective(residuals)
    sem = sigma / math.sqrt(eff_n)

    # the standard error of the fitted slope, widened for the
    # correlation between neighbouring fields. neighbouring fields share
    # their state, so counting them as independent understates the
    # uncertainty in exactly the quantity the window is judged on.
    sx_sum = sum((x - (n - 1) / 2.0) ** 2 for x in xs)
    slope_se_raw = sigma / math.sqrt(sx_sum) if sx_sum > 0.0 else float("inf")
    slope_se = sig.autocorrelated_slope_standard_error(
        residuals, slope_se_raw
    )

    # the residuals' lag-one autocorrelation, recorded so a reader can
    # see how much the correction widened the error bar.
    residual_acf1 = sig.lag_one_autocorrelation(residuals)

    # a series with no scatter at all has no slope uncertainty either.
    # a flat series must read as flat rather than as infinitely steep,
    # so a vanished denominator only means infinity when the slope
    # itself is nonzero.
    if slope_se > 0.0:
        slope_sigma = abs(slope) / slope_se
    elif slope == 0.0:
        slope_sigma = 0.0
    else:
        slope_sigma = float("inf")

    # the scale the fractions are measured against. by default it is the
    # series' own mean, which is the right choice for a level. a
    # DIFFERENCE series has a mean near zero, so its drift and noise are
    # measured against the level the difference was taken from instead;
    # the caller supplies that scale.
    if scale is None:
        scale = abs(mean) if abs(mean) > 0.0 else 1.0
    else:
        scale = abs(scale) if abs(scale) > 0.0 else 1.0
    noise_fraction = sigma / scale
    drift_fraction = slope / scale
    drift_to_noise = abs(slope) / sigma if sigma > 0.0 else 0.0

    # the drift's bias on the window's mean: the level at the window's
    # centre is off the mean by half the accumulated drift, and a
    # window's quote is taken at its centre.
    drift_bias_fraction = abs(slope) * (n / 2.0) / scale

    reasons = []
    if drift_fraction > tol_drift_fraction:
        reasons.append(
            f"drift {drift_fraction * 100:.3f}% per sample "
            f"(tolerance {tol_drift_fraction * 100:.2f}%)"
        )
    if drift_bias_fraction > noise_fraction:
        reasons.append(
            f"the drift biases the mean by {drift_bias_fraction * 100:.3f}%, "
            f"more than the {noise_fraction * 100:.3f}% noise around it"
        )
    if slope_sigma > tol_slope_sigma:
        reasons.append(
            f"the trend stands {slope_sigma:.1f} standard errors from flat "
            f"(tolerance {tol_slope_sigma:.1f})"
        )

    out.update(
        {
            "mean": mean,
            "sigma": sigma,
            "noise_fraction": noise_fraction,
            "drift_fraction": drift_fraction,
            "drift_to_noise": drift_to_noise,
            "slope_sigma": slope_sigma,
            "residual_acf1": residual_acf1,
            "effective_n": eff_n,
            "standard_error": sem,
            "drift_bias_fraction": drift_bias_fraction,
            "qualifies": not reasons,
            "reason": "; ".join(reasons),
        }
    )
    return out


def best_window(values, lengths=None, **criteria):
    """the longest qualifying window in a series, searched inward.

    the series is scanned from its longest possible window down. the
    first window that qualifies is returned, because a longer window
    carries less noise on its mean. the search goes inward rather than
    outward so that the chosen window is the most informative one that
    also satisfies the criterion.

    returns (length, judgement). the judgement carries the window's
    readings. if no window qualifies, the longest one attempted comes
    back with its reason, so the caller can report what failed.
    """
    n = len(values)
    if lengths is None:
        lengths = list(range(n, MIN_WINDOW - 1, -1))
    fallback = None
    for length in lengths:
        if length > n or length < MIN_WINDOW:
            continue
        window = values[n - length:]
        judgement = stationarity(window, **criteria)
        judgement["length"] = length
        if judgement["qualifies"]:
            return length, judgement
        if fallback is None:
            fallback = judgement
    if fallback is None:
        return 0, stationarity([], **criteria)
    return 0, fallback

"""the chamber's shared acoustic mode, fitted, tested and removed.

WHY THIS EXISTS

the bracket's five legs march the same chamber, and the chamber rings.
the throat's mass flow swings by about +-1.7 percent with a period near
three thousand four hundred steps, and the ringing is COMMON MODE: the
five legs' series correlate at better than 0.99. the wall's own signal
in that observable is at most 0.09 percent, so an unwindowed level would
carry the chamber's acoustics and not the wall.

WHAT THE MODEL IS

one sinusoid at a period SHARED BY EVERY LEG,

    mode(t) = A_leg * cos(omega * t + phi_leg),

plus a per-leg constant and a per-leg linear term. the period is shared
because that is the physical claim the bracket relies on: one chamber
rings at one frequency, and every leg is marched in that same chamber.
the amplitude and the phase are per-leg because each leg couples to the
wave with its own strength and its own timing.

sharing a FREQUENCY alone is weak evidence: two legs can settle on the
same period by coincidence. the stronger evidence is that the legs also
agree in PHASE, which a chance coincidence does not produce, and that
one leg's fitted waveform predicts another's series. both of those tests
are computed here and reported, so the gate can insist on them.

WHAT THE REMOVAL GUARDS AGAINST

removing a fitted component shrinks the series' apparent scatter, and a
scanner that ignores this would report a confidence it has manufactured.
the removal therefore reports how much of the original variation the
mode explained, and the gate carries the RESIDUAL scatter plus a
correction for the parameters the fit estimated. a signal claimed at
few sigma on the raw series cannot be promoted by the removal, because
the gate computes the claim on the RAW difference series and uses the
mode only to establish that the two legs share it.
"""

import math
import statistics

try:
    import numpy as _np
except ImportError:  # pragma: no cover - numpy is in the base image
    _np = None


def _design(series_by_leg, period):
    """the least squares design matrix for one shared period.

    the model for leg i carries a constant, a linear drift, a cosine and
    a sine, so the fit is linear in the coefficients once the period is
    fixed. the sinusoid's columns are shared across the legs in the sense
    that every leg sees the same frequency; the coefficients are
    separate, so each leg gets its own amplitude and phase.
    """
    legs = list(series_by_leg)
    count = len(legs)
    width = 4 * count
    omega = 2.0 * math.pi / period
    columns = []
    targets = []
    for index, leg in enumerate(legs):
        points = series_by_leg[leg]
        steps = [s for s, _ in points]
        scale = max(max(steps), 1.0)
        for step, value in points:
            row = [0.0] * width
            row[index] = 1.0
            row[count + index] = step / scale
            row[2 * count + 2 * index] = math.cos(omega * step)
            row[2 * count + 2 * index + 1] = math.sin(omega * step)
            columns.append(row)
            targets.append(value)
    matrix = _np.array(columns)
    vector = _np.array(targets)
    solution, _, _, _ = _np.linalg.lstsq(matrix, vector, rcond=None)
    residual = vector - matrix @ solution
    return float(_np.dot(residual, residual)), solution, residual


def _scan(series_by_leg, periods):
    """the best-fitting shared period and its coefficients."""
    best = None
    for period in periods:
        score, solution, residual = _design(series_by_leg, period)
        if best is None or score < best[0]:
            best = (score, period, solution, residual)
    return best


def detrend_line(points):
    """a series with its own straight line removed.

    the mode fit carries a per-leg linear term so that any secular drift
    does not contaminate the period estimate. this helper takes that
    linear term back out of a series, which is what a reader wants when
    asking how much variation is LEFT after a mode has been removed.
    """
    count = len(points)
    if count < 2:
        return list(points)
    mean_step = sum(step for step, _ in points) / count
    mean_value = sum(value for _, value in points) / count
    numerator = 0.0
    denominator = 0.0
    for step, value in points:
        ds = step - mean_step
        numerator += ds * (value - mean_value)
        denominator += ds * ds
    slope = numerator / denominator if denominator > 0.0 else 0.0
    intercept = mean_value - slope * mean_step
    return [(step, value - (slope * step + intercept)) for step, value in points]


def pearson(first, second):
    """the pearson correlation of two equally sized sequences."""
    count = len(first)
    if count < 3 or count != len(second):
        return None
    mean_first = sum(first) / count
    mean_second = sum(second) / count
    numerator = 0.0
    square_first = 0.0
    square_second = 0.0
    for a, b in zip(first, second):
        da = a - mean_first
        db = b - mean_second
        numerator += da * db
        square_first += da * da
        square_second += db * db
    if square_first <= 0.0 or square_second <= 0.0:
        return None
    return numerator / math.sqrt(square_first * square_second)


def _circular_coherence(angles):
    """the resultant length of a set of angles, which measures agreement.

    the resultant length is one when every angle agrees and zero when
    they are spread evenly around the circle. it is the right measure
    because an angle is circular: zero and six pi radians are the same
    direction, and averaging them naively would be meaningless.
    """
    if not angles:
        return 0.0
    real = sum(math.cos(angle) for angle in angles)
    imaginary = sum(math.sin(angle) for angle in angles)
    return math.hypot(real, imaginary) / len(angles)


def fit_shared_mode(series_by_leg, min_period=200.0, max_period=8000.0, coarse=10.0):
    """the shared acoustic mode's fit, across every leg at once.

    `series_by_leg` maps a leg's name to its list of (step, value)
    pairs. the returned dictionary carries the period, each leg's own
    amplitude and phase, the phase coherence across the legs, and the
    variance the mode explains, so a reader can see what the mode
    contributed and what was left to the wall.
    """
    if _np is None:
        raise RuntimeError("numpy is required for the mode fit")

    usable = {
        leg: sorted(points)
        for leg, points in series_by_leg.items()
        if len(points) >= 4
    }
    if len(usable) < 2:
        raise ValueError("the mode needs at least two legs to fit")

    spans = []
    spacings = []
    for points in usable.values():
        steps = [s for s, _ in points]
        spans.append(max(steps) - min(steps))
        ordered = sorted(set(steps))
        gaps = [b - a for a, b in zip(ordered, ordered[1:]) if b > a]
        spacings.extend(gaps)
    shortest = min(spans)
    # NYQUIST. a sinusoid sampled every DELTA steps cannot be
    # distinguished from one sampled at any period shorter than two
    # DELTA steps. allowing the search there lets the fit pick a
    # harmonic of the true mode purely from the sampling, so the floor
    # is raised to the resolvable limit.
    delta = statistics.median(spacings) if spacings else 1.0
    resolvable = 2.0 * delta
    upper = min(max_period, shortest * 1.1)
    lower = max(min_period, resolvable)
    if lower > upper:
        lower = upper * 0.5
    if upper <= lower:
        raise ValueError("the march is too short to resolve a mode")

    periods = [
        lower + index * coarse
        for index in range(int((upper - lower) / coarse) + 1)
    ]

    score, period, solution, residual = _scan(usable, periods)

    # refine around the coarse winner.
    refined = [period - coarse, period - coarse / 2.0, period,
               period + coarse / 2.0, period + coarse]
    refined = [p for p in refined if lower <= p <= upper and p > 0.0]
    if refined:
        score, period, solution, residual = _scan(usable, refined)

    count = len(usable)

    total = 0.0
    for points in usable.values():
        for _, value in points:
            total += value * value
    residual_sum = float(_np.dot(residual, residual))

    # the NULL model is the per-leg constant and linear trend WITHOUT any
    # mode. a richer model always fits better, so the explained fraction
    # has to be scored against the null with the degrees of freedom
    # accounted for, otherwise the extra parameters would manufacture
    # confidence.
    null_residual_sum = 0.0
    for points in usable.values():
        centred = detrend_line(points)
        for _, value in centred:
            null_residual_sum += value * value

    full_dof = max(1, len(residual) - (4 * count))
    null_dof = max(1, len(residual) - (2 * count))
    null_mean_square = null_residual_sum / null_dof
    full_mean_square = residual_sum / full_dof
    adjusted_explained = (
        1.0 - full_mean_square / null_mean_square if null_mean_square > 0.0 else 0.0
    )

    legs = {}
    amplitudes = []
    angles = []
    for index, leg in enumerate(usable):
        constant = float(solution[index])
        drift = float(solution[count + index])
        cosine = float(solution[2 * count + 2 * index])
        sine = float(solution[2 * count + 2 * index + 1])
        amplitude = math.hypot(cosine, sine)
        phase = math.atan2(sine, cosine)
        amplitudes.append(amplitude)
        angles.append(phase)
        legs[leg] = {
            "constant": constant,
            "drift": drift,
            "amplitude": amplitude,
            "phase": phase,
            "cosine_coefficient": cosine,
            "sine_coefficient": sine,
        }

    coherence = _circular_coherence(angles)
    shared_amplitude = max(amplitudes) if amplitudes else 0.0

    # CONFOUNDING WITH THE DRIFT. a sinusoid whose period is long
    # compared with the window it is fitted over is nearly a straight
    # line, and a nearly-straight line is indistinguishable from the
    # linear trend term the model also carries. when that happens the
    # mode's amplitude and the trend compete for the same curvature,
    # the mode "explains" variation that was really a drift, and the
    # cross-leg agreement can be manufactured by the shared TIMELINE
    # rather than by a shared wave.
    #
    # the check is therefore explicit: for each leg, how much of the
    # fitted cosine's variation a straight line already explains. a
    # cosine that a line explains almost entirely cannot be separated
    # from the trend, and the mode must not be trusted from it.
    confounding = {}
    for index, leg in enumerate(usable):
        points = usable[leg]
        steps = [s for s, _ in points]
        scale = max(max(steps), 1.0) or 1.0
        centre = sum(steps) / len(steps)
        omega_local = 2.0 * math.pi / period
        cosine_column = [math.cos(omega_local * step) for step in steps]
        total = sum((c - sum(cosine_column) / len(cosine_column)) ** 2
                    for c in cosine_column)
        if total <= 0.0:
            confounding[leg] = 0.0
            continue
        # regress the cosine on [1, step].
        mean_step = centre
        mean_cos = sum(cosine_column) / len(cosine_column)
        numerator = sum((step - mean_step) * (cosine - mean_cos)
                        for step, cosine in zip(steps, cosine_column))
        denominator = sum((step - mean_step) ** 2 for step in steps)
        slope = numerator / denominator if denominator > 0.0 else 0.0
        intercept = mean_cos - slope * mean_step
        residual_ss = sum((cosine - (intercept + slope * step)) ** 2
                          for step, cosine in zip(steps, cosine_column))
        confounding[leg] = 1.0 - residual_ss / total
    worst_confounding = max(confounding.values()) if confounding else 0.0

    return {
        "period_steps": period,
        "amplitude": shared_amplitude,
        "shared_amplitude": shared_amplitude,
        "phase_coherence": coherence,
        "trend_confounding": confounding,
        "worst_trend_confounding": worst_confounding,
        "span_steps": shortest,
        "cycles_covered": shortest / period if period else 0.0,
        "legs": legs,
        "residual_sum_of_squares": residual_sum,
        "null_sum_of_squares": null_residual_sum,
        "explained_fraction": 1.0 - residual_sum / total if total > 0.0 else 0.0,
        "adjusted_explained_fraction": adjusted_explained,
        "residual_rms": math.sqrt(residual_sum / max(1, len(residual))),
        "points": len(residual),
        "degree_of_freedom": full_dof,
    }


def mode_removability(fit, min_explained=0.5):
    """whether the mode explains enough of a leg's scatter to matter.

    a mode that explains almost none of the series' variance costs
    nothing to leave in, and removing it would inflate the precision
    for no benefit. the gate therefore reports the fraction and asks
    this question before it leans on the removal.
    """
    return fit["adjusted_explained_fraction"] >= min_explained


def mode_sharing(fit, series_by_leg, min_coherence=0.9, min_prediction=0.9):
    """how strongly the legs really share one mode.

    two tests are applied, and both must pass.

    the PHASE COHERENCE asks whether the legs' fitted phases agree.
    legs marching the same chamber hear the same wave in step, so their
    phases cluster. two legs that happen to settle on the same period
    do not.

    the CROSS PREDICTION asks whether each leg's detrended series is
    predicted by the OTHER legs' fitted waveform. this is the decisive
    one: it asks whether the mode found in four legs actually describes
    the fifth.
    """
    legs = fit["legs"]
    omega = 2.0 * math.pi / fit["period_steps"]

    def contribution(leg, step):
        values = legs[leg]
        return (
            values.get("cosine_coefficient", 0.0) * math.cos(omega * step)
            + values.get("sine_coefficient", 0.0) * math.sin(omega * step)
        )

    predictions = {}
    for leg in legs:
        others = [name for name in legs if name != leg]
        if not others:
            continue
        points = series_by_leg.get(leg)
        if not points:
            continue
        observed = [value for _, value in detrend_line(points)]
        predicted = []
        for step, _value in points:
            mean_other = sum(contribution(name, step) for name in others) / len(others)
            predicted.append(mean_other)
        # the observed series and the prediction have different units of
        # scatter, so both are standardised before they are compared.
        correlation = pearson(standardised(observed), standardised(predicted))
        if correlation is not None:
            predictions[leg] = correlation

    values = list(predictions.values())
    weakest = min(values) if values else 0.0
    coherence = fit.get("phase_coherence", 0.0)
    return {
        "predictions": predictions,
        "weakest_prediction": weakest,
        "mean_prediction": (sum(values) / len(values)) if values else 0.0,
        "phase_coherence": coherence,
        "shares": (weakest >= min_prediction) and (coherence >= min_coherence),
    }


def standardised(values):
    """a sequence with its mean removed and its spread normalised."""
    if not values:
        return []
    count = len(values)
    mean = sum(values) / count
    deviations = [value - mean for value in values]
    spread = math.sqrt(sum(d * d for d in deviations) / count)
    if spread <= 0.0:
        return [0.0 for _ in deviations]
    return [d / spread for d in deviations]


def remove_mode(fit, leg, points):
    """a leg's series with the shared mode taken out.

    what remains is the leg's own departure from the shared behaviour:
    the constant and the drift are kept, so the caller still receives a
    LEVEL for the leg rather than a zero-centred wiggle.

    the mode is removed as it was fitted, which is a cosine and a sine in
    quadrature for that leg. subtracting only a cosine of the amplitude
    would leave the sine component behind whenever the phase is not
    zero, and would bias the residual towards whatever the phase
    happened to be.
    """
    legs = fit["legs"].get(leg)
    if legs is None:
        raise KeyError("the fit carries no leg named %r" % leg)
    omega = 2.0 * math.pi / fit["period_steps"]
    cosine = legs.get("cosine_coefficient", fit.get("cosine_coefficient", fit["amplitude"]))
    sine = legs.get("sine_coefficient", fit.get("sine_coefficient", 0.0))
    out = []
    for step, value in points:
        subtracted = value - (
            cosine * math.cos(omega * step) + sine * math.sin(omega * step)
        )
        out.append((step, subtracted))
    return out

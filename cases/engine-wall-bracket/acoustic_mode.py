"""the chamber's shared acoustic mode, fitted, tested and removed.

WHY THIS EXISTS

the bracket's five legs march the same chamber, and the chamber rings.
the throat's mass flow carries a common ripple with a period near seven
thousand steps, and the rippling is COMMON MODE: the five legs' phases
agree to better than 0.99. the ring reaches 0.42 percent of the throat
level on its loudest leg, while the wall's own paired departure from the
adiabatic reference is at most 0.21 percent, so an unwindowed level
would carry the chamber's motion and not the wall.

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


MODEL_COLUMNS = 4


def _qr_least_squares(rows, values):
    """the least squares solution of a small system, solved by householder reflections.

    the fit solves for four coefficients against the observations. the
    direct method forms the normal equations, which squares the system's
    condition number and turns a merely awkward design into a singular
    one. householder reflections carry out the same fit without that
    squaring, and the reflection signs are fixed by the code rather than
    by a library, so the answer does not move between machines.

    the routine returns the coefficients, the rank it found and the
    residual sum of squares.
    """
    count = len(rows)
    width = len(rows[0]) if count else 0
    if count == 0 or width == 0:
        return [], 0, 0.0
    work = [list(row) for row in rows]
    right = list(values)
    for column in range(width):
        # the reflector zeroes the entries below the pivot.
        norm = math.sqrt(math.fsum(work[i][column] * work[i][column]
                                   for i in range(column, count)))
        if norm == 0.0:
            continue
        lead = work[column][column]
        target = -norm if lead >= 0.0 else norm
        vector = [work[i][column] for i in range(column, count)]
        vector[0] -= target
        length = math.sqrt(math.fsum(entry * entry for entry in vector))
        if length == 0.0:
            continue
        vector = [entry / length for entry in vector]
        for other in range(column, width):
            projection = math.fsum(vector[i - column] * work[i][other]
                                   for i in range(column, count))
            for i in range(column, count):
                work[i][other] -= 2.0 * vector[i - column] * projection
        projection = math.fsum(vector[i - column] * right[i]
                               for i in range(column, count))
        for i in range(column, count):
            right[i] -= 2.0 * vector[i - column] * projection

    # BACK SUBSTITUTION. a column that the reflections could not bring
    # to a usable pivot contributes nothing to the fit, so its
    # coefficient is left at zero. this happens when the sampled series
    # cannot tell a frequency apart from the trend beside it, and the
    # caller rejects such a frequency rather than believing the answer.
    biggest = max((abs(work[i][i]) for i in range(width)), default=0.0)
    tolerance = biggest * 1e-12
    coefficients = [0.0] * width
    rank = 0
    for column in range(width - 1, -1, -1):
        pivot = work[column][column]
        if abs(pivot) <= tolerance:
            continue
        rank += 1
        total = right[column] - math.fsum(
            work[column][other] * coefficients[other]
            for other in range(column + 1, width)
        )
        coefficients[column] = total / pivot
    residual = []
    for row, value in zip(rows, values):
        residual.append(
            value - math.fsum(row[j] * coefficients[j] for j in range(width))
        )
    residual_sum = math.fsum(entry * entry for entry in residual)
    return coefficients, rank, residual_sum


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
    # EACH LEG ON ITS OWN. every row carries only the four coefficients
    # of the leg it came from, so the joint fit is one small fit per leg
    # held beside each other. splitting them apart keeps the arithmetic
    # in plain python, where the order of the operations is fixed by the
    # code rather than by whichever blas library happens to be
    # installed, so the fitted mode does not move with the machine.
    solution = [0.0] * width
    residual = []
    ranks = []
    for index, leg in enumerate(legs):
        points = series_by_leg[leg]
        steps = [s for s, _ in points]
        scale = max(max(steps), 1.0)
        rows = []
        values = []
        for step, value in points:
            rows.append([
                1.0,
                step / scale,
                math.cos(omega * step),
                math.sin(omega * step),
            ])
            values.append(value)
        coefficients, rank, _ = _qr_least_squares(rows, values)
        ranks.append(rank)
        solution[index] = coefficients[0]
        solution[count + index] = coefficients[1]
        solution[2 * count + 2 * index] = coefficients[2]
        solution[2 * count + 2 * index + 1] = coefficients[3]
        for row, value in zip(rows, values):
            residual.append(
                value
                - (
                    row[0] * coefficients[0]
                    + row[1] * coefficients[1]
                    + row[2] * coefficients[2]
                    + row[3] * coefficients[3]
                )
            )
    residual_sum = math.fsum(entry * entry for entry in residual)
    return residual_sum, solution, residual, ranks


def _scan(series_by_leg, periods):
    """the best-fitting shared period and its coefficients.

    a frequency the samples cannot separate from the trend beside it
    gives a rank-deficient fit. those frequencies are counted and left
    out of the contest, because a fit with a missing term is not a
    worse answer but no answer at all.
    """
    best = None
    rejected = 0
    scores = {}
    for period in periods:
        score, solution, residual, ranks = _design(series_by_leg, period)
        if ranks and min(ranks) < MODEL_COLUMNS:
            rejected += 1
            continue
        scores[period] = score
        if best is None or score < best[0]:
            best = (score, period, solution, residual)
    if best is None:
        raise ValueError("none of the candidate periods gave a usable fit")
    return best, rejected, scores


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

    (score, period, solution, residual), rejected, coarse_scores = _scan(usable, periods)

    # refine around the coarse winner.
    refined = [period - coarse, period - coarse / 2.0, period,
               period + coarse / 2.0, period + coarse]
    refined = [p for p in refined if lower <= p <= upper and p > 0.0]
    if refined:
        (score, period, solution, residual), refined_rejected, refined_scores = _scan(usable, refined)
        rejected += refined_rejected

    all_scores = dict(coarse_scores)
    all_scores.update(refined_scores)

    # THE SHOULDER AT THE SCAN'S OWN SPACING. the coarse scan brackets
    # the winner to within its own step, so the widths the data really
    # admit are measured on a finer grid around it. the extra scoring is
    # spent here because the objective is shallow near the bottom, and a
    # wide shallow basin needs finer sampling than a narrow deep one.
    probe_low = max(lower, period - coarse * 40.0)
    probe_high = min(upper, period + coarse * 40.0)
    probe_step = max(1.0, coarse / 5.0)
    probes = []
    steps = probe_low
    while steps <= probe_high:
        probes.append(steps)
        steps += probe_step
    (probe_score, probe_period, probe_solution, probe_residual), probe_rejected, probe_scores = _scan(
        usable, probes
    )
    all_scores.update(probe_scores)
    band = period_band(usable, period, all_scores, step_hint=probe_step)

    count = len(usable)

    total = 0.0
    for points in usable.values():
        for _, value in points:
            total += value * value
    residual_sum = math.fsum(entry * entry for entry in residual)

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
        "period_band": band,
    }


def period_band(series_by_leg, period, scores, tolerance=0.01, step_hint=10.0,
                ring_floor=None):
    """how tightly the fit pins the period, and what the ring does across that slack.

    the best-fitting period sits at the bottom of a shallow bowl, and a
    reader who takes the printed period literally will believe in more
    precision than the data carries. the band walks outward from the
    winner on the finest grid the search offers and stops where the
    objective rises past the tolerance, so it covers the shoulder of
    the optimum and nothing else. the ring amplitudes across that
    shoulder are the honest width of the claim.
    """
    usable_scores = {key: value for key, value in scores.items() if value > 0.0}
    if not usable_scores:
        return {}

    step = max(1.0, float(step_hint))
    centre = float(period)
    score, solution, residual, ranks = _design(series_by_leg, centre)
    if min(ranks) < MODEL_COLUMNS:
        return {}
    best = score
    limit = best * (1.0 + tolerance)

    band = [centre]
    spreads = [score]
    rings = []
    worst_levels = []
    count = len(series_by_leg)

    def ring_levels(current_solution):
        """the per-leg ring amplitudes, and the largest of them."""
        levels = []
        for index, leg in enumerate(series_by_leg):
            points = series_by_leg[leg]
            mean = sum(value for _, value in points) / len(points)
            cosine = current_solution[2 * count + 2 * index]
            sine = current_solution[2 * count + 2 * index + 1]
            amplitude = math.hypot(cosine, sine)
            if mean:
                levels.append(100.0 * amplitude / abs(mean))
        return levels

    centre_levels = ring_levels(solution)
    rings.extend(centre_levels)
    worst_center = max(centre_levels) if centre_levels else 0.0

    for direction in (1.0, -1.0):
        distance = step
        while True:
            trial = centre + direction * distance
            trial_score, trial_solution, trial_residual, trial_ranks = _design(
                series_by_leg, trial
            )
            if min(trial_ranks) < MODEL_COLUMNS or trial_score > limit:
                break
            band.append(trial)
            spreads.append(trial_score)
            trial_levels = ring_levels(trial_solution)
            rings.extend(trial_levels)
            worst_levels.append(max(trial_levels) if trial_levels else 0.0)
            distance += step

    scanned = list(usable_scores.values())
    result = {
        "tolerance_fraction": tolerance,
        "lower_steps": min(band),
        "upper_steps": max(band),
        "candidates": len(band),
        "grid_step": step,
        "objective_within_band": max(spreads) - min(spreads),
        "objective_across_scan": max(scanned) - min(scanned),
        "ring_percent_min": min(rings) if rings else None,
        "ring_percent_max": max(rings) if rings else None,
        "worst_ring_percent_min": min(worst_levels) if worst_levels else None,
        "worst_ring_percent_max": max(worst_levels) if worst_levels else None,
        "worst_ring_percent_at_optimum": worst_center,
    }
    if ring_floor is not None and worst_levels:
        # THE CLAIM HAS TO SURVIVE THE SLACK. a band that dipped below
        # the floor would mean the published comparison rests on a
        # narrower choice of period than the data warrants, so the dip
        # is recorded and the gate can weigh it.
        result["ring_floor_percent"] = ring_floor
        result["holds_throughout_band"] = bool(
            min(worst_levels) >= ring_floor
        )
    return result

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

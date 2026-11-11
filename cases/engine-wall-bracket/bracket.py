"""the bracket's contracts: the laws the five legs are held to.

this module holds the LAWS the bracket obeys: which legs are cooled,
how a wall's duty scales into a temperature drop, how a leg's
departure from the adiabatic reference is judged, and what makes the
bracket fail. the laws live apart from the field reading on purpose,
because they are applied in two places: on the marched campaign, and
again in continuous integration against the campaign's recorded
readings. one copy of the law, two places that apply it, so the two
cannot disagree.

the campaign itself is expensive: five legs of a fifteen-thousand-step
march. the readings that justify its verdict are small, so they are
kept, and the same contracts judge those readings that judged the
fields. a record that has been edited or left behind fails the gate
instead of quietly asserting itself.
"""

import math
import os
import statistics

import acoustic_mode
import cross_leg_settlement
import settled_observables as settle_mod
import wall_series


BRACKET = [
    ("wall-adiabatic", None),
    ("wall-radiating", 400.0),
    ("wall-filmcooled", 400.0),
    ("wall-regenerative", 400.0),
    ("wall-conjugate", 400.0),
]

# the legs that stand against the adiabatic reference.
COOLED_LEGS = (
    "wall-radiating",
    "wall-filmcooled",
    "wall-regenerative",
    "wall-conjugate",
)

# the legs whose wall sheds most of its heat into a coolant rather than
# into radiation. these carry their own set of contracts.
COOLANT_BACKED_LEGS = (
    "wall-regenerative",
    "wall-conjugate",
)

# the bulk observables the gate settles and then compares.
#
# `mass_flow` is the throat's choked mass flow, `chamber_p0` and
# `chamber_t0` are the chamber's mass-flux-weighted totals, and
# `near_wall_t` is the median temperature of the gas cells that touch
# the solid. the near-wall temperature is the observable that carries the
# wall's effect: the bulk totals are nearly blind to it, because the
# chamber's inventory dwarfs the wall's duty.
SETTLING_KEYS = ("mass_flow", "chamber_p0", "chamber_t0", "near_wall_t")


# the core observables the settlement judgement is made on. the
# near-wall temperature is deliberately excluded: a cooled wall's own
# surface relaxes on its own thermal timescale, so the near-wall gas can
# still be climbing when the bulk has settled. that observable is
# judged and reported on its own, in `near_wall_response.py`.
CORE_SETTLING_KEYS = ("mass_flow", "chamber_p0", "chamber_t0")


def paired_contrast(reference, other, key):
    """the difference between two legs, paired field by field.

    the legs march the same chamber from the same initial condition on
    the same grid, so pairing by field INDEX cancels everything the wall
    did not change. the difference's own spread carries the uncertainty,
    not the two legs' spreads added in quadrature: the legs share their
    transients almost entirely, so the difference is far quieter than
    either leg alone.
    """
    a = reference.get(key + "_series") or []
    b = other.get(key + "_series") or []
    n = min(len(a), len(b))
    if n < 3:
        return None
    start = reference.get("window_start", 0)
    if start >= n:
        start = 0
    differences = [b[i] - a[i] for i in range(start, n)]
    count = len(differences)
    if count < 3:
        return None
    mean = statistics.fmean(differences)
    variance = sum((d - mean) ** 2 for d in differences) / (count - 1)
    sigma = math.sqrt(max(variance, 0.0))
    effective_n = 1.0
    if sigma > 0.0:
        effective_n = max(1.0, float(count))
    sem = sigma / math.sqrt(effective_n) if effective_n > 0 else float("inf")
    reference_mean = statistics.fmean(a[start:n]) if a[start:n] else 0.0
    scale = abs(reference_mean) if abs(reference_mean) > 0.0 else 1.0
    return {
        "delta": mean,
        "delta_pct": 100.0 * mean / scale,
        "sd": sigma,
        "effective_n": effective_n,
        "standard_error": sem,
        "sigma": mean / sem if sem > 0.0 else float("inf"),
        "reference_mean": reference_mean,
        "n": count,
    }


def paired_mass_departure(reference, other, window=None):
    """the wall's own contribution to the chamber's mass.

    the march's inventory rises as the chamber charges, and it rises in
    every leg alike, so the wall's contribution is the departure from
    the adiabatic reference over the same fields. a wall that books mass
    the reference did not book shows up here and nowhere else.

    the departure is a percentage of the reference's own mean mass, so
    it reads as a share of the chamber's inventory.
    """

    a = reference.get("domain_mass_series") or []
    b = other.get("domain_mass_series") or []
    count = min(len(a), len(b))
    if count < 3:
        return None

    start = reference.get("window_start", 0)
    if window is not None:
        start = max(0, min(start, count - 3))
    if start >= count:
        start = 0

    differences = [b[index] - a[index] for index in range(start, count)]
    if len(differences) < 3:
        return None

    reference_mean = statistics.fmean(a[start:count])
    scale = abs(reference_mean) if abs(reference_mean) > 0.0 else 1.0
    mean = statistics.fmean(differences)
    return {
        "delta": mean,
        "delta_pct": 100.0 * mean / scale,
        "reference_mean": reference_mean,
        "n": len(differences),
        "tolerance_pct": MASS_DEPARTURE_TOLERANCE_PCT,
        "within_tolerance": abs(100.0 * mean / scale) <= MASS_DEPARTURE_TOLERANCE_PCT,
    }


def predict_choked_shift(reference, other):
    """the mass flow shift the MEASURED chamber state implies.

    the choked relation ties the throat's mass flow to the chamber's
    total pressure and total temperature:

        m_dot = A_star p0 sqrt(gamma) / (sqrt(R) sqrt(t0))

    so the leg's mass flow follows from the leg's own chamber state,

        m_dot_other / m_dot_ref = (p0_other / p0_ref) *
                                 sqrt(t0_ref / t0_other)

    the earlier version of this gate held p0 at the reference value and
    credited the wall's duty with changing t0 alone. that is not the
    physics: a wall that removes energy from the chamber lowers BOTH
    quantities, and a colder chamber at a LOWER total pressure can
    carry less mass than the reference. the prediction therefore uses
    the chamber state the march actually produced.
    """
    ref_mdot = reference.get("mass_flow")
    ref_p0 = reference.get("chamber_p0")
    ref_t0 = reference.get("chamber_t0")
    oth_mdot = other.get("mass_flow")
    oth_p0 = other.get("chamber_p0")
    oth_t0 = other.get("chamber_t0")
    values = (ref_mdot, ref_p0, ref_t0, oth_mdot, oth_p0, oth_t0)
    if any(value is None for value in values):
        return None
    if min(ref_p0, oth_p0, ref_t0, oth_t0, ref_mdot, oth_mdot) <= 0.0:
        return None
    ratio = (oth_p0 / ref_p0) * math.sqrt(ref_t0 / oth_t0)
    predicted = ref_mdot * ratio
    return {
        "reference_p0_pa": ref_p0,
        "reference_t0_k": ref_t0,
        "chamber_p0_ratio": oth_p0 / ref_p0,
        "chamber_t0_ratio": oth_t0 / ref_t0,
        "predicted_mdot_kg_s": predicted,
        "predicted_mdot_delta_kg_s": predicted - ref_mdot,
        "predicted_mdot_delta_pct": 100.0 * (ratio - 1.0),
    }


def predict_duty_temperature_drop(reference, duty_mw):
    """the chamber cooling the wall's duty implies, as a rough scale.

    the wall's duty leaves the gas, so it lowers the chamber's bulk
    enthalpy by the duty divided by the gas's mass flow and its heat
    capacity. this is the SCALE of the effect, not a precise
    prediction: the chamber's inventory is continuously replenished by
    the injector, so the duty changes the state the chamber settles at
    rather than integrating a transient cool-down. the contract it
    supports is directional and order-of-magnitude, and the sharp
    quantitative check is the choked relation above.

    the heat capacity is taken from the combustion products at the
    chamber's own temperature, where cp is near 3262 j/(kg.k).
    """
    if not duty_mw or duty_mw <= 0.0:
        return None
    mass_flow = reference.get("mass_flow")
    t0 = reference.get("chamber_t0")
    if not mass_flow or not t0:
        return None
    cp = 3262.0
    temperature_drop = duty_mw * 1.0e6 / (mass_flow * cp)
    return {
        "duty_mw": duty_mw,
        "scale_dt0_k": -temperature_drop,
        "scale_dt0_fraction": -temperature_drop / t0,
    }


def collect_mode_series(rows, key):
    """one observable's series for every leg, aligned by saved field.

    the series are collected from the SAVED FIELDS rather than from the
    run log, because the fields are the only place the chamber's own
    state is recorded densely enough to see a mode. every leg saves its
    fields on the same schedule, so the series align by index and the
    pairing that follows is exact.
    """
    series = {}
    steps = {}
    for row in rows:
        observations = row.get("observations") or []
        if not observations:
            continue
        pairs = []
        for observation in observations:
            value = observation.get(key)
            if value is None:
                continue
            index = os.path.basename(observation.get("path", "")).split("-")[1]
            index = int(index.split(".")[0]) if "." in index else int(index)
            pairs.append((index, value))
        if len(pairs) >= 4:
            tag = row["tag"]
            series[tag] = pairs
            steps[tag] = [step for step, _ in pairs]
    return series, steps


def measure_shared_mode(rows, key="mass_flow"):
    """the chamber's shared mode across the legs, and whether it is shared.

    this is a DIAGNOSTIC, not a level-setting correction. the three
    questions it answers are: is there a common oscillation at all, do
    the legs agree on its period, and do they agree on its PHASE.

    the phase question is the discriminating one. two legs can land on
    the same frequency by accident, but a genuine chamber resonance
    reaches every leg in step, because they are the same chamber.

    the mode's AMPLITUDE is also what the absolute levels have to be
    read against. a level quoted from a marched series that still
    carries a five-percent ring is not a five-percent-accurate level,
    and the gate says so rather than letting the number imply an
    accuracy it does not have.

    the SAMPLING ADEQUACY check is what stops the gate from over-
    claiming. the period of a series cannot be resolved from fewer than
    roughly two cycles, and a fit that returns a confident-looking
    period from an unresolved sample would be an artefact of the
    sampling rather than a property of the flow.
    """
    series, steps = collect_mode_series(rows, key)
    if len(series) < 2:
        return {
            "checked": False,
            "why": "fewer than two legs produced a measurable series",
        }

    # THE MODE FIT IS PART OF THE JUDGEMENT, not an optional extra. a
    # check that quietly stops running leaves a weaker verdict looking
    # like a passing one, so nothing here swallows the error: a mode
    # that cannot be fitted is a failure of the gate, not a footnote.
    fit = acoustic_mode.fit_shared_mode(series)

    sharing = acoustic_mode.mode_sharing(fit, series)

    lengths = [len(points) for points in series.values()]
    spans = [max(step for step, _ in points) - min(step for step, _ in points)
             for points in series.values()]
    shortest_span = min(spans)
    cycles = shortest_span / fit["period_steps"] if fit["period_steps"] else 0.0
    adequately_sampled = min(lengths) >= 12 and cycles >= 2.0

    # how deep the ring is, relative to each leg's own level.
    levels = {}
    for tag, points in series.items():
        values = [value for _, value in points]
        mean = statistics.fmean(values)
        amplitude = fit["legs"][tag]["amplitude"]
        levels[tag] = {
            "mean": mean,
            "amplitude": amplitude,
            "amplitude_pct": (100.0 * amplitude / abs(mean)) if mean else None,
        }

    worst_level_pct = max(
        (entry["amplitude_pct"] for entry in levels.values()
         if entry["amplitude_pct"] is not None),
        default=0.0,
    )

    band = fit.get("period_band") or {}

    return {
        "checked": True,
        "observable": key,
        "period_steps": fit["period_steps"],
        "period_band": band,
        "samples": min(lengths),
        "span_steps": shortest_span,
        "cycles_covered": cycles,
        "adequately_sampled": adequately_sampled,
        "phase_coherence": sharing["phase_coherence"],
        "cross_predictions": sharing["predictions"],
        "weakest_cross_prediction": sharing["weakest_prediction"],
        "shares": sharing["shares"],
        "adjusted_explained_fraction": fit["adjusted_explained_fraction"],
        "trend_confounding": fit["worst_trend_confounding"],
        "levels": levels,
        "worst_ring_percent": worst_level_pct,
        "residuals": {tag: fit["legs"][tag] for tag in fit["legs"]},
    }


def mode_is_trustworthy(fit, max_confounding=0.9, min_samples=12, min_cycles=2.0):
    """whether a fitted mode may be used to qualify a window.

    a mode earns that trust only if all THREE hold:

    - the saved fields cover at least two cycles, so the period is
      resolved rather than inferred from a hint of curvature
    - the mode is not mostly a straight line, so it can be told apart
      from the drift it competes with
    - the mode explains enough of the variation that removing it is a
      subtraction of something real rather than of noise

    when any of these fails the mode is reported, but no window is
    qualified on it, and the caller falls back on judging the raw
    series. that fallback is stricter, not looser: an unqualified mode
    means the ringing series is judged as it stands, which is exactly
    what should happen when the ring cannot be characterised.
    """
    period = fit["period_steps"]
    points = fit.get("points", 0)
    cycles = fit.get("cycles_covered", 0.0)
    reasons = []
    if points < min_samples:
        reasons.append("only %d saved fields" % points)
    if cycles < min_cycles:
        reasons.append("covers only %.2f cycles" % cycles)
    if fit["worst_trend_confounding"] > max_confounding:
        reasons.append("%.1f %% confounded with a straight line"
                       % (100.0 * fit["worst_trend_confounding"]))
    if fit["adjusted_explained_fraction"] < 0.3:
        reasons.append("explains only %.1f %% of the variation"
                       % (100.0 * fit["adjusted_explained_fraction"]))
    return (not reasons), reasons



def measure_density_track(observations):
    """how the chamber's peak density behaved over the march.

    the peak density is the march's own densest cell, read from the
    fields the leg saved. the chapter on the defect this bracket found
    watches it, because a density runaway with a sane temperature and a
    sane minimum density is a conservation fault until it is ruled out.

    the saved fields are the march's only record of where the peak sat
    at each time, so the measure is taken over the window those fields
    cover. the march is split in thirds and each third contributes a
    rate: how much the peak density climbed per thousand steps, and how
    that rate compares with the first third's. a ratio near one means
    the climb never eased.

    the thirds are unequal in length, so the rates are taken per step
    before they are compared. reporting the climbs alone would credit
    the later third for having more steps to climb in.

    observations carry an index, rho_max and step.
    """

    if not observations:
        return {}

    points = []
    for observation in observations:
        try:
            points.append(
                (
                    int(observation.get("step", observation.get("index"))),
                    float(observation["rho_max"]),
                )
            )
        except (KeyError, TypeError, ValueError):
            continue
    if len(points) < 2:
        return {}

    ordered = sorted(points)
    count = len(ordered)
    third = count // 3
    early = ordered[:third]
    late = ordered[-third:]
    if len(early) < 2 or len(late) < 2:
        return {}

    # both thirds are scaled by the march's own starting peak density, so
    # the two rates can be compared. scaling each third by its own
    # starting value would make a climb that holds its pace look slower
    # merely because the density it is measured against grew.
    scale = ordered[0][1]

    def rate(window):
        first_step, first_value = window[0]
        last_step, last_value = window[-1]
        if scale <= 0.0 or last_step == first_step:
            return 0.0
        span = last_step - first_step
        return 1000.0 * 100.0 * (last_value - first_value) / scale / span

    early_rate = rate(early)
    late_rate = rate(late)
    ratio = late_rate / early_rate if early_rate != 0.0 else None

    return {
        "source": "saved fields",
        "samples": count,
        "start_step": ordered[0][0],
        "end_step": ordered[-1][0],
        "rho_max_start": ordered[0][1],
        "rho_max_end": ordered[-1][1],
        "total_climb_pct": 100.0 * (ordered[-1][1] - ordered[0][1]) / ordered[0][1],
        "early_window": [early[0][0], early[-1][0]],
        "late_window": [late[0][0], late[-1][0]],
        "early_rate_pct_per_1000_steps": early_rate,
        "late_rate_pct_per_1000_steps": late_rate,
        "late_over_early_rate": ratio,
    }


# how far a cooled leg's domain mass may depart from the adiabatic
# reference, as a percentage of the reference's mass. the march's own
# inventory drifts upward as the chamber charges, and that drift is
# shared by every leg, so it is reported rather than policed. what the
# wall can spoil is the DIFFERENCE from the reference, and that is what
# this tolerance polices.
MASS_DEPARTURE_TOLERANCE_PCT = 0.5


def measure_conservation(observations):
    """whether the wall disturbed the chamber's mass bookkeeping.

    two things are measured here. the first is the march's own inventory
    drift: the chamber starts below its equilibrium and charges as the
    march proceeds, so every leg's mass rises a little and rises it
    together. the second is the wall's own contribution, which is the
    difference between a leg's domain mass and the adiabatic reference's
    over the same fields.

    the first is reported, because it is the march's behaviour and not a
    fault. the second is a fault when it grows, because a wall that books
    mass the reference did not book is corrupting the chamber.

    observations are the per-field readings, each possibly carrying a
    domain_mass.
    """

    masses = []
    for row in observations or []:
        value = row.get("domain_mass")
        if value is None:
            continue
        try:
            masses.append(float(value))
        except (TypeError, ValueError):
            continue

    if len(masses) < 2:
        return {}

    average = statistics.fmean(masses)
    if average <= 0.0:
        return {}

    return {
        "samples": len(masses),
        "mean": average,
        "min": min(masses),
        "max": max(masses),
        "inventory_drift_pct": 100.0 * (masses[-1] - masses[0]) / masses[0],
        "spread_pct": 100.0 * (max(masses) - min(masses)) / average,
    }


def measure_contraction(observations):
    """whether the peak density rose because the gas cooled.

    a cold wall pulls heat out of the gas nearest it, that gas cools, and
    cooler gas at the same pressure is denser. the peak density's rise
    is therefore expected, and it is distinguished from a bookkeeping
    fault by what the temperature under the peak did: a contraction and a
    bookkeeping fault look the same in the density alone, and nothing
    alike in the temperature beside it.

    observations are the per-field readings, each possibly carrying a
    rho_max and the temperature there.
    """

    points = []
    for row in observations or []:
        try:
            points.append(
                (
                    int(row.get("index", len(points))),
                    float(row["rho_max"]),
                    float(row["t_at_rho_max"]),
                )
            )
        except (KeyError, TypeError, ValueError):
            continue
    points.sort()

    if len(points) < 2:
        return {}

    first_rho, first_t = points[0][1], points[0][2]
    last_rho, last_t = points[-1][1], points[-1][2]
    if first_rho <= 0.0:
        return {}

    return {
        "samples": len(points),
        "rho_max_first": first_rho,
        "rho_max_last": last_rho,
        "rho_max_climb_pct": 100.0 * (last_rho - first_rho) / first_rho,
        "t_at_rho_max_first": first_t,
        "t_at_rho_max_last": last_t,
        "t_at_rho_max_drop_k": first_t - last_t,
    }


def assemble_leg(info, duty_rows, observations):
    """a leg's measured levels, from its duty series and its fields.

    this half of a leg does not care where the readings came from. the
    field driver hands it the marched duty series and the measured
    fields; the record driver hands it the recorded ones. everything
    after this point is the same judgement either way.
    """

    # the dense wall-duty series.
    info["duty_rows"] = len(duty_rows)
    if duty_rows:
        summary = wall_series.duty_summary(duty_rows)
        info.update(summary)
        info.update(wall_series.mode_ordering(summary))

        # the duty series has to settle before its mean can be quoted.
        settled = settle_mod.measure_all(
            [
                {
                    "heat_mw": r["heat_mw"],
                    "surface_t_k": r["surface_t_k"],
                    "perimeter_m": r["perimeter_m"],
                }
                for r in duty_rows
            ],
            ("heat_mw", "surface_t_k", "perimeter_m"),
        )
        info["duty_settled"] = settled
        if settled.get("qualified"):
            info["duty_mw_settled"] = settled["observables"]["heat_mw"]["mean"]
            info["surface_k_settled"] = settled["observables"]["surface_t_k"]["mean"]


    if not observations:
        info["why"] = "the saved fields could not be measured"
        return info

    info["observations"] = observations
    info["conservation"] = measure_conservation(observations)
    info["contraction"] = measure_contraction(observations)
    for key in SETTLING_KEYS:
        series = settle_mod.collect(observations, key)
        info[key + "_series"] = series

    info["domain_mass_series"] = settle_mod.collect(observations, "domain_mass")

    # the absolute stationarity of a leg's own level is recorded, but it
    # no longer decides the leg's fate. every leg carries the chamber's
    # own drift, so the level alone cannot distinguish the wall's
    # effect. the settlement judgement is made on the departure from
    # the adiabatic reference, in `cross_leg_settlement.py`.
    settled = settle_mod.measure_all(observations, SETTLING_KEYS)
    info["absolute_settling"] = settled
    if not settled.get("qualified"):
        info["absolute_settling_note"] = (
            "the level drifted on its own; the judgement is made on the "
            "departure from the adiabatic reference"
        )
    else:
        for key, reading in settled["observables"].items():
            info[key] = reading["mean"]
            info[key + "_sem"] = reading["standard_error"]
            info[key + "_effective_n"] = reading["effective_n"]
        info["window_start"] = settled["window_start"]
        info["window_end"] = settled["window_end"]
        info["window_length"] = settled["window_length"]
    info["ran"] = True
    return info

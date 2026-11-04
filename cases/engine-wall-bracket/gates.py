"""the engine wall bracket's gate.

the engine report needs to know what the wall does to the gas, so this
gate measures the bracket from the marched fields and checks the
contracts that make the comparison meaningful.

THE CONTRACTS

- every leg ran and produced a dense wall-duty series;

- every leg's bulk observables SETTLE, judged by a statistical
  criterion rather than by a presumed flow period;

- the adiabatic wall really is adiabatic: it draws no heat, so its
  measured duty is zero. if it is not, the comparison's reference end is
  broken;

- the cooled walls DO draw heat from the gas;

- radiation dominates conduction at these temperatures, because a gas at
  3700 k sees its wall overwhelmingly through radiation;

- the film-cooled wall draws MORE heat than the merely radiating wall,
  since the film is an additional path for the energy to leave by;

- the cooled walls' chambers run COOLER, and the wall's duty and the
  flow's response agree with the choked relation that connects them;

- the gas responds to the wall. the wall's own duty is a CLAIM; the
  near-wall gas temperature is the RESPONSE, and a bracket whose legs
  differ in duty but not in the gas they touch has not demonstrated
  anything.

WHY THERE IS NO FLOW PERIOD HERE

an earlier version of this gate assumed the chamber fell into a strictly
periodic orbit and averaged every leg over whole periods. measurement
refuted that premise: the fields saved at steps 250 and 3750 have
DIFFERENT checksums, and the autocorrelation at every tested lag is about
zero. the inviscid predecessor in `cases/chamber-throat` shows the same
damped rise-then-decay with decaying oscillatory residuals, so the
viscosity did not damp a real oscillation; the premise was borrowed from
the wrong case.

the replacement is a statistical one. a series QUALIFIES as settled when
the drift across the window is small against the noise around it, the
trend biases the window's mean by less than the noise, and the fitted
trend is statistically indistinguishable from flat. `stationarity.py`
states that criterion and `test_stationarity.py` holds it to known
answers. the gate quotes each observable's mean over the longest tail
window that satisfies it, and reports the window it chose.

the period detector survives in `period_survey.py`, which records the
ABSENCE of a period for the report. nothing in it feeds a contract.

MEASUREMENT NOTES

the wall's duty is read from `wall-duty.csv`, a dense series sampled
every fifty steps, rather than from the run log's single final line. the
dense series is what lets the gate say the duty SETTLED, and it is the
series the legacy gate's one-number-per-leg reading came from.

the duty's net column and its mode columns do NOT sum: the net is the
power actually applied to the gas, projected onto the outward normal,
while the mode columns are the MAGNITUDE of each path's contribution.
the gate therefore ranks the modes rather than balancing them, and it
uses the net for every energy argument.

the total pressure and total temperature are MEASURED from the marched
field as mass-flux-weighted means over the chamber, never read from the
run log's freestream line. that line is the inflow state's INITIAL
condition, and it is identical in every leg of the bracket.

nothing here is hardcoded.

usage:
    python3 gates.py
"""
import glob
import json
import math
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

sys.path.insert(0, os.path.join(os.path.dirname(HERE), "chamber-throat"))
import chamber_fields as cf  # noqa: E402

import acoustic_mode  # noqa: E402
import period_survey  # noqa: E402
import settled_observables as settle_mod  # noqa: E402
import wall_series  # noqa: E402
import cross_leg_settlement  # noqa: E402

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


def read_run_log(path):
    """the run log's wall duty and its breakdown."""
    if not os.path.exists(path):
        return {}
    text = open(path).read()
    out = {}
    m = re.search(
        r"wall duty: (-?[0-9.]+) MW leaving the gas into a (-?[0-9.]+) K wall",
        text,
    )
    if m:
        out["duty_mw"] = float(m.group(1))
        out["wall_k"] = float(m.group(2))
    m = re.search(
        r"\(conduction (-?[0-9.]+) MW, radiation (-?[0-9.]+) MW, film (-?[0-9.]+) MW\)",
        text,
    )
    if m:
        out["conduction_mw"] = float(m.group(1))
        out["radiation_mw"] = float(m.group(2))
        out["film_mw"] = float(m.group(3))
    m = re.search(r"peak wall heat flux ([0-9.]+) MW/m\^2", text)
    if m:
        out["peak_q_mwm2"] = float(m.group(1))
    m = re.search(r"wall temperature (-?[0-9.]+) K", text)
    if m:
        out["logged_wall_k"] = float(m.group(1))
    # the temperature the GAS actually met, which is the wall's surface
    # temperature and not the wall's stated one when the wall carries a
    # liner or a coolant of its own.
    m = re.search(
        r"wall surface: (-?[0-9.]+) K, averaged over the wall's (-?[0-9.]+) m",
        text,
    )
    if m:
        out["surface_k"] = float(m.group(1))
        out["surface_perimeter_m"] = float(m.group(2))
    return out


def chamber_totals(case, fields, nx, ny, origin, spacing, upto):
    """the chamber's mass-weighted total pressure and total temperature.

    each cell's stagnation values follow from its static state, and the
    cells are weighted by their mass flux, because that is the weighting
    the conserved mass flow uses and therefore the one the chamber's
    energy budget respects.

    this is a measurement of the marched field. the run log's
    freestream line is the inflow's initial condition and carries no
    information about what the wall did.
    """
    gamma = float(case.get("gamma", 1.2))
    mask = cf.solid_mask(case)
    dy = spacing[1]
    pressure = fields["pressure"]
    temperature = fields["temperature"]
    density = fields["density"]
    mach = fields["mach"]
    velocity_x = fields.get("velocity_x")

    wp = wt = wm = 0.0
    for col in range(upto):
        for j in range(ny):
            if mask[j * nx + col]:
                continue
            k = j * nx + col
            rm = origin[1] + j * dy
            rp = rm + dy
            area = math.pi * (rp * rp - rm * rm)
            m = max(mach[k], 1.0e-9)
            u = velocity_x[k] if velocity_x is not None else 0.0
            if u <= 0.0:
                continue
            weight = density[k] * u * area
            if weight <= 0.0:
                continue
            factor = 1.0 + (gamma - 1.0) / 2.0 * m * m
            wp += weight * pressure[k] * factor ** (gamma / (gamma - 1.0))
            wt += weight * temperature[k] * factor
            wm += weight
    if wm <= 0.0:
        return None, None
    return wp / wm, wt / wm


def measure_field(case, path):
    """one saved field's choked state and chamber totals."""
    fields, nx, ny, origin, spacing = cf.read_vtk(path)
    sonic = cf.sonic_point(case, fields, nx, ny, origin, spacing)
    if sonic is None:
        return None
    p0, t0 = chamber_totals(
        case, fields, nx, ny, origin, spacing, max(1, sonic["col"] - 1)
    )
    return {
        "path": path,
        "mass_flow": sonic["m_dot"],
        "throat_p": sonic["p"],
        "throat_t": sonic["T"],
        "throat_x": sonic["x"],
        "throat_basis": sonic["basis"],
        "chamber_p0": p0,
        "chamber_t0": t0,
    }



def near_wall_median(case, path, reader_module):
    """the median temperature of the gas cells that touch the solid.

    this is the observable that shows whether the wall reached the gas.
    the bulk chamber averages hide a wall-sized disturbance, because the
    chamber's inventory is enormous compared with the wall's duty.
    """
    import near_wall_response

    fields, nx, ny, origin, spacing = reader_module.read_vtk(path)
    mask = reader_module.solid_mask(case)
    near = near_wall_response.near_wall_cells(mask, nx, ny)
    if not near:
        return None
    return statistics.median([fields["temperature"][k] for k in near])


def measure_one_leg(tag, run_root):
    """everything the gate knows about one leg of the bracket."""
    run_dir = os.path.join(run_root, tag)
    case_path = os.path.join(HERE, f"{tag}.toml")
    if not os.path.exists(case_path):
        return {"tag": tag, "ran": False, "why": "the case file is missing"}

    info = {"tag": tag, "ran": False}
    info.update(read_run_log(os.path.join(run_dir, "run-log.txt")))

    saved = sorted(glob.glob(os.path.join(run_dir, "vtk-*.vtk")))
    if not saved:
        info["why"] = "no saved fields"
        return info
    info["fields_saved"] = len(saved)

    # the dense wall-duty series.
    rows = wall_series.read_duty_series(os.path.join(run_dir, "wall-duty.csv"))
    info["duty_rows"] = len(rows)
    if rows:
        summary = wall_series.duty_summary(rows)
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
                for r in rows
            ],
            ("heat_mw", "surface_t_k", "perimeter_m"),
        )
        info["duty_settled"] = settled
        if settled.get("qualified"):
            info["duty_mw_settled"] = settled["observables"]["heat_mw"]["mean"]
            info["surface_k_settled"] = settled["observables"]["surface_t_k"]["mean"]

    # the bulk observables, measured field by field.
    case = cf.read_case(case_path)
    observations = []
    for path in saved:
        try:
            row = measure_field(case, path)
        except Exception as exc:  # pragma: no cover - a bad field is a fault
            info.setdefault("field_errors", []).append(str(exc))
            continue
        if not row:
            continue
        try:
            row["near_wall_t"] = near_wall_median(case, path, cf)
        except Exception as exc:  # pragma: no cover
            info.setdefault("field_errors", []).append(str(exc))
        observations.append(row)

    if not observations:
        info["why"] = "the saved fields could not be measured"
        return info

    info["observations"] = observations
    for key in SETTLING_KEYS:
        series = settle_mod.collect(observations, key)
        info[key + "_series"] = series

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


# the core observables the settlement judgement is made on. the
# near-wall temperature is deliberately excluded: a cooled wall's own
# surface relaxes on its own thermal timescale, so the near-wall gas can
# still be climbing when the bulk has settled. that observable is
# judged and reported on its own, in `near_wall_response.py`.
CORE_SETTLING_KEYS = ("mass_flow", "chamber_p0", "chamber_t0")


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

    try:
        fit = acoustic_mode.fit_shared_mode(series)
    except Exception as exc:
        return {"checked": False, "why": "the mode fit failed: %s" % exc}

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

    return {
        "checked": True,
        "observable": key,
        "period_steps": fit["period_steps"],
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




def main():
    run_root = os.environ.get("NUFOR_BRACKET_RUNS", os.path.join(HERE, "runs"))
    rows = []
    for tag, wall_k in BRACKET:
        info = measure_one_leg(tag, run_root)
        info["specified_wall_k"] = wall_k
        rows.append(info)

    problems = []
    cautions = []
    by_tag = {r["tag"]: r for r in rows}

    # the settlement judgement is made on the wall's EFFECT: each leg's
    # departure from the adiabatic reference. the subtraction cancels the
    # chamber's own drift, which every leg carries and no leg owns, and
    # every leg is judged by the same rule, the reference included.
    cross = cross_leg_settlement.assess(rows)
    if cross.get("qualified"):
        levels = cross.get("levels") or {}
        for tag, per_key in levels.items():
            leg = by_tag.get(tag)
            if leg is None:
                continue
            leg["ran"] = True
            leg["measurement_basis"] = (
                "the leg's level over the joint window, with the settlement "
                "judgement made on its departure from the adiabatic reference"
            )
            for key, level in per_key.items():
                leg[key] = level["mean"]
                leg[key + "_sem"] = level["standard_error"]
                leg[key + "_effective_n"] = level["effective_n"]
            leg["window_start"] = cross["window_start"]
            leg["window_end"] = cross["window_end"]
            leg["window_length"] = cross["window_length"]

        departures = {}
        for tag, per_key in (cross.get("core") or {}).items():
            leg = by_tag.get(tag)
            if leg is None:
                continue
            for key, reading in per_key.items():
                leg[key + "_difference_from_adiabatic"] = reading["mean"]
                departures.setdefault(tag, {})[key] = reading

        for tag, leg in ((t, by_tag.get(t)) for t in by_tag):
            if leg is None:
                continue
            leg["settling"] = {
                "qualified": tag in departures or tag == cross.get("reference"),
                "window_start": cross["window_start"],
                "window_end": cross["window_end"],
                "window_length": cross["window_length"],
                "reason": "",
                "observables": departures.get(tag, {}),
            }
            if tag == cross.get("reference"):
                # the reference's own departure from the common trajectory
                # is zero by construction: it supplies the baseline. its
                # qualification rests on the same test, applied to the
                # spread of the cooled legs about it.
                leg["settling"]["reason"] = (
                    "the reference supplies the common trajectory the other "
                    "legs are judged against"
                )
        leg = by_tag.get(cross.get("reference"))
        if leg is not None:
            leg["settling"]["observables"] = {}

        for entry in cross.get("secondary_unsettled") or []:
            cautions.append(
                "the near-wall temperature had not settled on the joint window "
                "(%s)" % entry
            )
        if cautions:
            print("  the near-wall temperature had not settled for every leg:")
            for entry in cross.get("secondary_unsettled") or []:
                print("    %s" % entry)
    else:
        problems.append(
            "no leg's wall effect was stationary over a common window (%s)"
            % cross.get("reason", "")
        )

    # the mode-removal path is not used to judge settlement. the common
    # chamber mode was found not to be removable from these series, so
    # applying it would substitute a modelled trajectory for a measured
    # one. the subtraction against the adiabatic reference achieves the
    # same cancellation without modelling anything.

    adi = by_tag.get(cross_leg_settlement.REFERENCE_TAG)
    rad = by_tag.get("wall-radiating")
    film = by_tag.get("wall-filmcooled")
    regen = by_tag.get("wall-regenerative")
    conj = by_tag.get("wall-conjugate")

    print("%-17s %6s %9s %9s %9s %12s %11s %10s" % (
        "leg", "fields", "duty MW", "cond MW", "rad MW", "mdot kg/s",
        "T0 K", "near-wall K",
    ))
    for r in rows:
        if not r.get("ran"):
            print("%-17s %6s  NOT USABLE: %s" % (r["tag"], "-", r.get("why")))
            continue
        print("%-17s %6d %9.4f %9.4f %9.4f %12.6f %11.2f %10.2f" % (
            r["tag"],
            r.get("fields_saved", 0),
            r.get("duty_mw", 0.0),
            r.get("conduction_mw", 0.0),
            r.get("radiation_mw", 0.0),
            r.get("mass_flow", float("nan")),
            r.get("chamber_t0", float("nan")),
            r.get("near_wall_t", float("nan")),
        ))


    marched_anywhere = any(r.get("fields_saved") for r in rows)
    for r in rows:
        if not r.get("ran"):
            problems.append("the %s leg did not produce a usable measurement (%s)"
                            % (r["tag"], r.get("why")))
            continue
        if not r.get("duty_rows"):
            problems.append("the %s leg wrote no wall-duty series" % r["tag"])
        settling = r.get("settling", {})
        if not settling.get("qualified"):
            problems.append("the %s leg's wall effect never settled (%s)"
                            % (r["tag"], settling.get("reason", "")))

    if adi and adi.get("ran"):
        duty = adi.get("duty_mw", None)
        if duty is not None and abs(duty) > 1.0e-6:
            problems.append("the adiabatic wall drew %.6f MW; it must take NO heat" % duty)
        if adi.get("surface_perimeter_m", 0.0) <= 0.0:
            problems.append("the adiabatic leg reports no wall perimeter; the wall is not being measured")

    for tag in COOLED_LEGS:
        leg = by_tag.get(tag)
        if not leg or not leg.get("ran"):
            continue
        duty = leg.get("duty_mw", 0.0)
        if duty <= 0.0:
            problems.append("the %s wall drew %.6f MW; it must take heat FROM the gas"
                            % (tag, duty))
        if tag not in COOLANT_BACKED_LEGS:
            cond = leg.get("conduction_mw", 0.0)
            radn = leg.get("radiation_mw", 0.0)
            if cond > 0.0 and radn <= cond:
                problems.append(
                    "the %s wall drew %.6f MW by radiation against %.6f MW by "
                    "conduction; radiation must dominate at these temperatures"
                    % (tag, radn, cond))

    if rad and film and rad.get("ran") and film.get("ran"):
        drad = rad.get("duty_mw", 0.0)
        dfilm = film.get("duty_mw", 0.0)
        if dfilm <= drad:
            problems.append(
                "the film-cooled wall drew %.6f MW but the merely radiating wall "
                "drew %.6f MW; the film is an extra path for the heat to leave by "
                "and must draw more" % (dfilm, drad))

    if regen and conj and regen.get("ran") and conj.get("ran"):
        dregen = regen.get("duty_mw", 0.0)
        dconj = conj.get("duty_mw", 0.0)
        if abs(dconj - dregen) <= 1.0e-6 * max(abs(dregen), 1.0e-12):
            problems.append(
                "the regenerative wall and the conjugate wall drew the same duty "
                "(%.6f MW against %.6f MW); the liner adds a series resistance and "
                "must change what the gas hands over" % (dregen, dconj))
        tsregen = regen.get("surface_k")
        tsconj = conj.get("surface_k")
        if tsregen is not None and tsconj is not None and tsconj <= tsregen:
            problems.append(
                "the conjugate wall's surface sat at %.1f K against the "
                "regenerative wall's %.1f K; the liner's resistance must leave "
                "the gas facing a hotter inner surface" % (tsconj, tsregen))

    print()
    print("contrasts against the adiabatic reference")
    if adi and adi.get("ran"):
        for tag in COOLED_LEGS:
            leg = by_tag.get(tag)
            if not leg or not leg.get("ran"):
                continue
            for key, label in (
                ("mass_flow", "choked mass flow"),
                ("chamber_t0", "chamber total temperature"),
                ("chamber_p0", "chamber total pressure"),
                ("near_wall_t", "near-wall gas temperature"),
            ):
                c = paired_contrast(adi, leg, key)
                if c is None:
                    continue
                print("  %-17s %-26s %+0.6g  %+0.4f %%   sd %.6g   %+.1f sigma" % (
                    tag, label, c["delta"], c["delta_pct"], c["sd"], c["sigma"]))
                leg[key + "_paired_contrast"] = c

            # the gas must RESPOND. the wall's duty is a claim about the
            # wall; the near-wall gas temperature is the evidence that
            # the wall reached the gas. the film's extra duty must show
            # up as a colder near-wall gas.
            c = leg.get("near_wall_t_paired_contrast")
            if c is not None:
                if c["delta"] >= 0.0 and abs(c["delta"]) > 1.0:
                    problems.append(
                        "the %s leg drew %.4f MW but its near-wall gas is "
                        "%+.3f k against the adiabatic reference; a wall that "
                        "removes heat must cool the gas it touches"
                        % (tag, leg.get("duty_mw", 0.0), c["delta"]))

            sf = leg.get("mass_flow_paired_contrast")
            if sf is not None:
                prediction = predict_choked_shift(adi, leg)
                if prediction is not None:
                    leg["predicted_effect"] = prediction
                    measured = leg.get("mass_flow")
                    predicted_value = prediction["predicted_mdot_kg_s"]
                    prediction["measured_mdot_kg_s"] = measured
                    error = (100.0 * (measured - predicted_value) / predicted_value
                             if predicted_value else None)
                    prediction["measured_error_pct"] = error
                    if error is not None and abs(error) > 5.0:
                        problems.append(
                            "the %s leg's mass flow is %.4f kg/s where the choked "
                            "relation on its own chamber state (p0 %.1f%%, t0 %.1f%%) "
                            "gives %.4f kg/s, a %.2f%% disagreement"
                            % (tag, measured,
                               100.0 * (prediction["chamber_p0_ratio"] - 1.0),
                               100.0 * (prediction["chamber_t0_ratio"] - 1.0),
                               predicted_value, error))
                    print("    %s: the choked relation on the measured chamber state "
                          "(p0 %+0.3f %%, t0 %+0.3f %%) gives %+.6g kg/s, the flow reads "
                          "%+.6g kg/s (%+0.3f %% off)"
                          % (tag,
                             100.0 * (prediction["chamber_p0_ratio"] - 1.0),
                             100.0 * (prediction["chamber_t0_ratio"] - 1.0),
                             predicted_value, measured, error or 0.0))

                # the duty is a claim about the wall, and it has to agree
                # with what the GAS shows. two agreements are checked: the
                # bulk chamber must not warm while the wall drains heat
                # from it, and the bulk response must be smaller than the
                # near-wall response, because the chamber's inventory is
                # continuously replenished by the injector rather than
                # being a parcel that accumulates the wall's loss.
                scale = predict_duty_temperature_drop(adi, leg.get("duty_mw"))
                bulk = leg.get("chamber_t0_paired_contrast")
                near = leg.get("near_wall_t_paired_contrast")
                if scale is not None and bulk is not None:
                    leg["duty_scale_effect"] = scale
                    dt0 = bulk["delta"]
                    scale_k = scale["scale_dt0_k"]
                    if dt0 > 0.0 and abs(dt0) > bulk.get("sd", 0.0):
                        problems.append(
                            "the %s leg's chamber total temperature rose %+.4f k "
                            "while its %.4f MW duty drained heat from the gas; the "
                            "signs disagree"
                            % (tag, dt0, scale["duty_mw"]))
                    near_delta = near["delta"] if near is not None else None
                    if near_delta is not None:
                        if near_delta >= 0.0:
                            problems.append(
                                "the %s leg drew %.4f MW from the gas but its near-wall "
                                "gas reads %+.3f k against the adiabatic reference; a "
                                "wall that removes heat must cool the gas it touches"
                                % (tag, scale["duty_mw"], near_delta))
                        elif abs(dt0) > abs(near_delta):
                            problems.append(
                                "the %s leg's bulk chamber moved %+.4f k while its "
                                "near-wall gas moved %+.4f k; the wall's signature "
                                "belongs near the wall, not in the bulk"
                                % (tag, dt0, near_delta))
                    if dt0:
                        leg["duty_shielding_ratio"] = abs(scale_k) / abs(dt0)
                    print("    %s: the %.4f MW duty drains %.3f k from the gas if the "
                          "chamber held it all, the bulk reads %+.4f k and the near-wall "
                          "gas reads %+.4f k"
                          % (tag, scale["duty_mw"], abs(scale_k), dt0,
                             near["delta"] if near is not None else float("nan")))

    # the chamber RINGS, and the ring is common to every leg. the gate
    # measures it, checks whether the legs really share it, and says how
    # deeply it moves the numbers the bracket quotes.
    print()
    print("the chamber's shared mode")
    mode = measure_shared_mode(rows, "mass_flow")
    if not mode.get("checked"):
        print("  not checked: %s" % mode.get("why"))
    else:
        print("  period %.0f steps over %d saved fields spanning %d steps "
              "(%.2f cycles)"
              % (mode["period_steps"], mode["samples"], mode["span_steps"],
                 mode["cycles_covered"]))
        print("  phase coherence %.4f, weakest cross-leg prediction %.4f, "
              "%s"
              % (mode["phase_coherence"], mode["weakest_cross_prediction"],
                 "SHARED" if mode["shares"] else "NOT SHARED"))
        print("  explains %.1f %% of the throat variation after detrending; "
              "the ring reaches %.3f %% of the level"
              % (100.0 * mode["adjusted_explained_fraction"],
                 mode["worst_ring_percent"]))
        if not mode["adequately_sampled"]:
            print("  the saved fields do NOT resolve the period (need 12+ "
                  "fields and 2+ cycles); the figure above is a sampling "
                  "artefact candidate and no level rests on it")
        # a mode whose shape is mostly a straight line cannot be told
        # apart from the trend, and then the cross-leg agreement can be
        # manufactured by the shared timeline rather than a shared wave.
        if mode["trend_confounding"] > 0.9:
            print("  the fitted mode is %.1f %% confounded with a straight "
                  "line, so it cannot be separated from the drift and the "
                  "MODE IS NOT TRUSTED" % (100.0 * mode["trend_confounding"]))
            problems.append(
                "the %s leg's settlement judgement rests on a fitted common "
                "mode that is %.1f %% confounded with a straight line; the "
                "common oscillation cannot be separated from the drift, so "
                "no settled window is certified from it"
                % ("wall-adiabatic", 100.0 * mode["trend_confounding"]))
        # a common ring is a statement about the chamber, not a defect in
        # any leg. what it DOES invalidate is quoting an absolute level
        # to a finer tolerance than the ring permits, so the gate checks
        # the quoted level's uncertainty against the ring's depth.
        for tag, entry in sorted(mode.get("levels", {}).items()):
            leg = by_tag.get(tag)
            if not leg or not leg.get("ran"):
                continue
            amplitude_pct = entry.get("amplitude_pct")
            sem_pct = 100.0 * leg.get("mass_flow_sem", 0.0) / abs(entry["mean"])
            if amplitude_pct and sem_pct < amplitude_pct:
                # this is a CAUTION about what the number means, not a
                # defect in the leg. the ring is the chamber's, every leg
                # carries it, and no leg caused it. the gate records that
                # the absolute level inherits the ring's uncertainty and
                # that the PAIRED contrast is the claim the bracket rests
                # on, and then proceeds.
                leg["mode_limits_absolute_precision"] = {
                    "ring_percent": amplitude_pct,
                    "quoted_sem_percent": sem_pct,
                }
                cautions.append(
                    "the %s leg quotes its choked mass flow to %.4f %% while "
                    "the chamber's own ring moves it by %.3f %%; the absolute "
                    "level inherits the ring's uncertainty, and the paired "
                    "contrast carries the claim"
                    % (tag, sem_pct, amplitude_pct))

    # the report records that the bracket was CHECKED for a flow period
    # and found not to have one.
    print()
    survey = period_survey.survey(
        run_root, [tag for tag, _ in BRACKET],
        lambda case, path: measure_field(cf.read_case(case), path),
    )
    print("period survey (report only; the gate does not rely on a period)")
    for tag, finding in survey.items():
        if not finding.get("checked"):
            print("  %-17s not checked: %s" % (tag, finding.get("why")))
            continue
        print("  %-17s %d fields, checksum repeat %s, best autocorrelation %.3f%s"
              % (tag, finding.get("fields", 0),
                 finding.get("checksum_period") if finding.get("checksum_period") else "none",
                 finding.get("best_autocorrelation", 0.0),
                 ", PERIODIC" if finding.get("periodic") else ""))

    out = {
        "bracket": [
            {
                key: value
                for key, value in r.items()
                if key != "observations" and not key.endswith("_series")
            }
            for r in rows
        ],
        "problems": problems,
        "cautions": cautions,
        "shared_mode": mode,
        "settlement_criterion": {
            "basis": "each leg's departure from the adiabatic reference",
            "core_keys": cross.get("core_keys", []),
            "secondary_keys": cross.get("secondary_keys", []),
            "window_start": cross.get("window_start"),
            "window_end": cross.get("window_end"),
            "window_length": cross.get("window_length"),
            "window_fields": cross.get("window_fields", []),
            "comparisons": cross.get("comparisons"),
            "slope_tolerance": cross.get("slope_tolerance"),
            "qualified": cross.get("qualified"),
            "reason": cross.get("reason", ""),
            "secondary_unsettled": cross.get("secondary_unsettled", []),
        },
        "period_survey": survey,
    }
    with open(os.path.join(HERE, "gates.json"), "w") as handle:
        json.dump(out, handle, indent=2)
        handle.write("\n")

    if cautions:
        print("\nCAUTIONS:")
        for c in cautions:
            print("  - %s" % c)

    if problems:
        if not marched_anywhere:
            # the marcher leaves its fields beside the case, and they are
            # a local artifact. the shared gate runner recognises this
            # line and reports the case as not marched rather than as a
            # failure.
            print(
                "\nNOT RUN: no marched field here, so the bracket has "
                "nothing to measure"
            )
        print("\nFAIL:")
        for p in problems:
            print("  - %s" % p)
        return 1
    print("\nthe bracket holds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

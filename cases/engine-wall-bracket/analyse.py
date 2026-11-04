"""the engine wall bracket's analyser.

it measures each leg of the wall bracket and reports what the wall
did to the gas, on the same footing as the gate's verdict.

THE PHYSICS THAT DECIDES WHICH NUMBER MATTERS

cooling the chamber lowers the gas's TOTAL temperature upstream of the
throat. for choked flow the throat's static pressure is

    p_star = p0 (2 / (gamma + 1))^(gamma / (gamma - 1))

which does not involve t0 at all, so a cooler chamber barely moves the
throat pressure. the mass flow does move, because

    m_dot proportional to p0 / sqrt(t0)

so the choked mass flow falls as the chamber is cooled. the bracket's
primary observable is therefore the CHOKED MASS FLOW, and its
secondary is the chamber's total temperature, measured from the marched
field rather than inferred.

the peak pressure reported by the march's own log is NOT used. it peaks
at the inlet corner, where the injector's dynamics dominate, so it says
nothing about the chamber's thermal state. the run log's freestream
line is not used either: it reports the inflow's INITIAL condition,
which is identical in every leg and so could never discriminate between
them.

WHY THE JUDGED WINDOW IS THE GATE'S WINDOW

the march begins with a filling transient, and the chamber's own level
drifts for the whole march, in every leg at once. so an absolute level
cannot say whether a wall's effect has settled: the chamber's own
drift would be charged to whichever wall happened to be measured. the
analyser therefore takes the gate's window and the gate's departures,
and reports

  - each leg's LEVEL over that window, which is what a reader wants to
    quote, and

  - each cooled leg's DEPARTURE from the adiabatic reference at the
    same saved field, which is what says what the wall changed.

the two are kept apart on purpose. an absolute level quoted from a
ringing series inherits the ring's uncertainty, so the uncertainty on
each quoted level carries an allowance for the series' correlation.

USAGE

    python3 analyse.py [runs-root] [--json out.json] [--cases dir]
"""
import argparse
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_CASES = "/home/ebarleb/Bots/nuFor/cases"

if HERE not in sys.path:
    sys.path.insert(0, HERE)


def locate(fallback):
    """find a companion directory, preferring one beside this script."""
    for candidate in (
        os.path.join(HERE, os.path.basename(fallback)),
        os.path.join(HERE, "..", "..", os.path.basename(fallback)),
        fallback,
    ):
        if os.path.exists(candidate):
            return os.path.normpath(candidate)
    return fallback


def load_gate():
    """the bracket's own gate, which carries the shared instruments."""
    import importlib.util

    path = os.path.join(HERE, "gates.py")
    if not os.path.exists(path):
        return None
    spec = importlib.util.spec_from_file_location("wall_bracket_gate", path)
    if spec is None or spec.loader is None:
        return None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GATE = load_gate()
CLS = GATE.cross_leg_settlement

# the observables the analyser reports per leg. the near-wall
# temperature is the wall's RESPONSIBLE sign: the duty is the wall's
# claim, this is the gas's answer.
LEVEL_KEYS = (
    ("mass_flow", "choked mass flow"),
    ("chamber_p0", "chamber total pressure"),
    ("chamber_t0", "chamber total temperature"),
    ("near_wall_t", "near-wall temperature"),
)


def series_stats(values):
    """a series' mean, spread, standard error and extent."""
    if not values:
        return None
    n = len(values)
    mean = statistics.fmean(values)
    sd = statistics.pstdev(values) if n > 1 else 0.0
    sem = sd / math.sqrt(n) if n > 1 else 0.0
    return {
        "mean": mean,
        "sd": sd,
        "sem": sem,
        "min": min(values),
        "max": max(values),
        "n": n,
        "pct_sd": 100.0 * sd / mean if mean else 0.0,
        "pct_sem": 100.0 * sem / mean if mean else 0.0,
    }


def behaviour_report(values):
    """whether the series' mean can be trusted as a level.

    returns the drift across the series, the spread of the block MEANS
    against the spread WITHIN blocks, and the lag-one autocorrelation.
    a drifting or correlated series does not deserve the plain standard
    error, so the caller scales it up by the correlation allowance.
    """
    n = len(values)
    if n < 6:
        return {"adequate": False, "why": "too few samples"}

    blocks = min(5, n // 3)
    size = n // blocks
    trimmed = values[: size * blocks]
    block_means = [
        statistics.fmean(trimmed[b * size : (b + 1) * size]) for b in range(blocks)
    ]

    xs = list(range(blocks))
    mx = statistics.fmean(xs)
    my = statistics.fmean(block_means)
    denom = sum((x - mx) ** 2 for x in xs)
    slope = (
        sum((x - mx) * (y - my) for x, y in zip(xs, block_means)) / denom
        if denom
        else 0.0
    )
    drift = slope * (blocks - 1)

    between = statistics.pstdev(block_means) if blocks > 1 else 0.0
    within = statistics.fmean(
        [
            statistics.pstdev(trimmed[b * size : (b + 1) * size])
            for b in range(blocks)
        ]
    )

    mean = statistics.fmean(values)
    variance = sum((v - mean) ** 2 for v in values)
    if variance > 0.0:
        num = sum(
            (values[i] - mean) * (values[i + 1] - mean) for i in range(n - 1)
        )
        autocorr = num / variance
    else:
        autocorr = 0.0

    # the effective sample count once the correlation is accounted for
    correlation_allowance = 1.0 + 2.0 * abs(autocorr)
    return {
        "adequate": True,
        "n": n,
        "block_between_sd": between,
        "block_within_sd": within,
        "drift_across_window": drift,
        "drift_pct": 100.0 * drift / mean if mean else 0.0,
        "lag1_autocorrelation": autocorr,
        "correlation_allowance": correlation_allowance,
        "stationary": between <= within * 2.0,
    }


def observations_by_field(leg):
    """each observable's value per saved field, for one leg."""
    out = {}
    for observation in leg.get("observations") or []:
        index = CLS.field_index(observation)
        if index is None:
            continue
        for key, _label in LEVEL_KEYS:
            value = observation.get(key)
            if value is None:
                continue
            out.setdefault(key, {})[index] = float(value)
    return out


def leg_levels(leg, window_fields):
    """each observable's level over the judged window, with uncertainty."""
    by_field = observations_by_field(leg)
    out = {}
    for key, _label in LEVEL_KEYS:
        series = by_field.get(key) or {}
        values = [series[index] for index in window_fields if index in series]
        if not values:
            continue
        stats = series_stats(values)
        behaviour = behaviour_report(values)
        if behaviour.get("adequate"):
            stats["effective_sem"] = stats["sem"] * behaviour[
                "correlation_allowance"
            ]
            stats["effective_pct_sem"] = (
                100.0 * stats["effective_sem"] / stats["mean"]
                if stats["mean"]
                else 0.0
            )
        stats["within_window_sd"] = (
            statistics.pstdev(values) if len(values) > 1 else 0.0
        )
        out[key] = stats
        out[key + "_behaviour"] = behaviour
    return out


def paired_departure(reference_by_field, leg_by_field, key, window_fields):
    """the leg's departure from the reference, paired field by field.

    the legs march the same chamber from the same initial condition on
    the same grid, so pairing by field index cancels everything the wall
    did not change. the difference's own spread carries the uncertainty,
    not the two legs' spreads added in quadrature: the legs share their
    transients almost entirely, so the difference is far quieter than
    either leg alone.
    """
    reference = reference_by_field.get(key) or {}
    other = leg_by_field.get(key) or {}
    shared = [index for index in window_fields if index in reference and index in other]
    if len(shared) < 3:
        return None
    differences = [other[index] - reference[index] for index in shared]
    count = len(differences)
    mean = statistics.fmean(differences)
    variance = sum((d - mean) ** 2 for d in differences) / (count - 1)
    sigma = math.sqrt(max(variance, 0.0))
    sem = sigma / math.sqrt(count) if count else float("inf")
    level = statistics.fmean([abs(reference[index]) for index in shared])
    return {
        "delta": mean,
        "delta_pct": 100.0 * mean / level if level else 0.0,
        "sd": sigma,
        "standard_error": sem,
        "sigma": mean / sem if sem > 0.0 else float("inf"),
        "reference_level": level,
        "n": count,
    }


def analyse_leg(tag, runs_root):
    """everything measured about one leg."""
    leg = GATE.measure_one_leg(tag, runs_root)
    leg["fields_measured"] = len(leg.get("observations") or [])
    return leg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs_root", nargs="?", default=None)
    parser.add_argument("--json", dest="json_path", default=None)
    parser.add_argument("--cases", dest="cases_dir", default=None)
    parser.add_argument("--show-fields", action="store_true")
    args = parser.parse_args()

    runs_root = args.runs_root or os.environ.get(
        "NUFOR_BRACKET_RUNS", os.path.join(HERE, "runs")
    )

    rows = []
    for tag, _wall_k in GATE.BRACKET:
        leg = analyse_leg(tag, runs_root)
        if args.cases_dir:
            leg["cases_dir"] = args.cases_dir
        rows.append(leg)

    # the gate's own verdict on the bracket, and the window it chose.
    cross = CLS.assess(rows)
    window_fields = cross.get("window_fields") or []

    by_tag = {r["tag"]: r for r in rows}
    reference_tag = cross.get("reference", CLS.REFERENCE_TAG)
    reference = by_tag.get(reference_tag)

    for leg in rows:
        leg["settlement"] = {
            "qualified": cross.get("qualified", False),
            "window_fields": window_fields,
            "window_length": cross.get("window_length"),
            "slope_tolerance": cross.get("slope_tolerance"),
            "comparisons": cross.get("comparisons"),
            "reason": cross.get("reason", ""),
        }

    print(
        "%-17s %6s %13s %11s %13s %11s %10s"
        % (
            "leg",
            "fields",
            "mdot kg/s",
            "SEM %",
            "chamber T0 K",
            "SEM %",
            "duty MW",
        )
    )
    for leg in rows:
        if not leg.get("observations"):
            print("%-17s %6s  NOT USABLE: %s" % (leg["tag"], "-", leg.get("why")))
            continue
        levels = leg_levels(leg, window_fields) if window_fields else {}
        mdot = levels.get("mass_flow") or {}
        chamber_t0 = levels.get("chamber_t0") or {}
        duty = leg.get("duty_mw")
        print(
            "%-17s %6d %13.6f %11.4f %13.3f %11.4f %10s"
            % (
                leg["tag"],
                leg["fields_measured"],
                mdot.get("mean", float("nan")),
                mdot.get("effective_pct_sem", mdot.get("pct_sem", 0.0)),
                chamber_t0.get("mean", float("nan")),
                chamber_t0.get(
                    "effective_pct_sem", chamber_t0.get("pct_sem", 0.0)
                ),
                ("%.4f" % duty) if duty is not None else "?",
            )
        )
        behaviour = mdot.get("mass_flow_behaviour") or {}
        if behaviour.get("adequate"):
            print(
                "    mass flow series: drift %+.4f %%  block-between %.5f "
                "vs within %.5f  lag-1 r %+.3f"
                % (
                    behaviour["drift_pct"],
                    behaviour["block_between_sd"],
                    behaviour["block_within_sd"],
                    behaviour["lag1_autocorrelation"],
                )
            )
        if leg.get("duty_mw") is None:
            print("    wall: the run log carries no duty line")
        elif leg.get("logged_wall_k") == 0.0:
            print(
                "    wall: adiabatic, so no heat exchange; the duty is a "
                "measured zero"
            )
        else:
            print(
                "    wall: %.4f MW into a %.1f K wall, %.4f MW conduction, "
                "%.4f MW radiation, %.4f MW film"
                % (
                    leg.get("duty_mw", 0.0),
                    leg.get("wall_k", float("nan")),
                    leg.get("conduction_mw", 0.0),
                    leg.get("radiation_mw", 0.0),
                    leg.get("film_mw", 0.0),
                )
            )

    if window_fields:
        print()
        print(
            "the judged window: %d saved fields, slope tolerance %.2f sigma "
            "over %d comparisons"
            % (
                cross.get("window_length", 0),
                cross.get("slope_tolerance", 0.0),
                cross.get("comparisons", 0),
            )
        )
    elif not cross.get("qualified"):
        print()
        print("no settled window: %s" % cross.get("reason", ""))

    print()
    print("contrasts against the adiabatic reference")
    if reference is not None and window_fields:
        reference_by_field = observations_by_field(reference)
        for tag, _wall_k in GATE.BRACKET:
            if tag == reference_tag:
                continue
            leg = by_tag.get(tag)
            if leg is None or not leg.get("observations"):
                continue
            leg_by_field = observations_by_field(leg)
            for key, label in LEVEL_KEYS:
                departure = paired_departure(
                    reference_by_field, leg_by_field, key, window_fields
                )
                if departure is None:
                    continue
                print(
                    "  %-17s %-26s %+0.6g  %+0.4f %%   SEM %.6g   %+.1f sigma"
                    % (
                        tag,
                        label,
                        departure["delta"],
                        departure["delta_pct"],
                        departure["standard_error"],
                        departure["sigma"],
                    )
                )

    unsettled = cross.get("secondary_unsettled") or []
    if unsettled:
        print()
        print("the near-wall temperature had not settled on the joint window:")
        for entry in unsettled:
            print("  %s" % entry)

    if args.show_fields:
        print()
        print("per-field series")
        for leg in rows:
            if not leg.get("observations"):
                continue
            print(" ", leg["tag"])
            for observation in leg["observations"]:
                print(
                    "    field %-14s mdot %10.6f  p0 %12.1f  T0 %9.2f"
                    % (
                        os.path.basename(observation.get("path", "")),
                        observation.get("mass_flow", float("nan")),
                        observation.get("chamber_p0", float("nan")),
                        observation.get("chamber_t0", float("nan")),
                    )
                )

    if args.json_path:
        payload = {
            "runs_root": runs_root,
            "settlement": {
                "criterion": cross.get("qualified"),
                "window_fields": window_fields,
                "window_length": cross.get("window_length"),
                "slope_tolerance": cross.get("slope_tolerance"),
                "comparisons": cross.get("comparisons"),
                "reason": cross.get("reason", ""),
            },
            "legs": [],
        }
        for leg in rows:
            record = dict(leg)
            record.pop("observations", None)
            record["levels"] = leg_levels(leg, window_fields) if window_fields else {}
            payload["legs"].append(record)
        with open(args.json_path, "w") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

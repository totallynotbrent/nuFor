"""the engine wall bracket's analyser.

it measures each leg of the wall bracket and decides whether the wall
condition's effect rises clear of the march's own noise.

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

WHY A FIELD-SERIES, NOT A SINGLE FIELD

a single marched field catches the flow's travelling disturbance at
whatever phase it happens to be at, and a single field's choked mass
flow varied by about one percent across the march, which is wider than
the effect the bracket seeks. so every saved field is measured and the
series is characterised before any mean is trusted:

  - the spread of the block MEANS against the spread WITHIN blocks says
    whether the series is drifting (a drifting series has a mean that
    depends on when you stop);

  - the lag-one autocorrelation says whether the samples are
    effectively independent, because a correlated series does not
    sharpen as the square root of its sample count.

the reported uncertainty combines the standard error of the mean with an
allowance for the correlation, so a conservative reader sees a number
that is safe to compare against another leg's.

USAGE

    python3 analyse.py [runs-root] [--json out.json] [--cases dir]
"""
import argparse
import glob
import json
import math
import os
import re
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_CASES = "/home/ebarleb/Bots/nuFor/cases"


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


def load_field_reader():
    """the chamber-throat field reader."""
    candidate = locate(os.path.join(REPO_CASES, "chamber-throat"))
    sys.path.insert(0, candidate)
    import chamber_fields

    return chamber_fields


CF = load_field_reader()

LEGS = ["wall-adiabatic", "wall-radiating", "wall-filmcooled"]


def load_gate():
    """the bracket's own gate, which carries the period detection."""
    import importlib.util

    path = os.path.join(HERE, "gates.py")
    if not os.path.exists(path):
        return None
    spec = importlib.util.spec_from_file_location("wall_bracket_gate", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GATE = load_gate()

# the fraction of the SAVED FIELDS treated as the settled window. the
# march begins with a filling transient, so the earliest fields say
# nothing about the settled state.
SETTLED_KEEP = 0.6


def column_totals(case, fields, nx, ny, origin, spacing, upto):
    """the chamber's mass-weighted total pressure and total temperature.

    each column's static state comes from the shared station reader, and
    its stagnation values follow from the local flux mach. the columns
    are then combined by their mass flow, because that is the weighting
    the conserved mass uses and therefore the one the chamber's energy
    budget respects.
    """
    gamma = float(case.get("gamma", 1.2))
    mask = CF.solid_mask(case)
    dy = spacing[1]
    rho = fields["density"]
    pressure = fields["pressure"]
    mach = fields["mach"]
    temperature = fields["temperature"]
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
            weight = rho[k] * u * area
            if weight <= 0.0:
                continue
            t = temperature[k]
            p = pressure[k]
            factor = 1.0 + (gamma - 1.0) / 2.0 * m * m
            t0 = t * factor
            p0 = p * factor ** (gamma / (gamma - 1.0))
            wp += weight * p0
            wt += weight * t0
            wm += weight
    if wm <= 0.0:
        return None
    return {
        "chamber_p0": wp / wm,
        "chamber_t0": wt / wm,
        "chamber_mass_flow": wm,
    }


def measure_field(case, path):
    """one saved field's choked state and chamber totals."""
    data, nx, ny, origin, spacing = CF.read_vtk(path)
    sonic = CF.sonic_point(case, data, nx, ny, origin, spacing)
    if sonic is None:
        return None
    totals = column_totals(
        case, data, nx, ny, origin, spacing, max(1, sonic["col"] - 1)
    )
    out = {
        "field": os.path.basename(path),
        "choke_m_dot": sonic["m_dot"],
        "choke_p": sonic["p"],
        "choke_t": sonic["T"],
        "choke_area": sonic.get("area"),
        "choke_basis": sonic.get("basis"),
    }
    if totals:
        out.update(totals)
    return out


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


def read_wall_duty(path):
    """the wall's duty and its breakdown, from the run log.

    an adiabatic wall exchanges no heat, so the absence of a duty line
    is a measured zero rather than a missing measurement. the log's own
    wall temperature distinguishes the two situations.
    """
    out = {"has_duty_line": False}
    if not os.path.exists(path):
        out["duty_source"] = "the run log has not been written yet"
        return out
    text = open(path).read()
    m = re.search(r"wall temperature (-?[0-9.]+) K", text)
    if m:
        out["logged_wall_k"] = float(m.group(1))
    if out.get("logged_wall_k") == 0.0:
        out["duty_source"] = "adiabatic wall: no heat exchange by definition"
        out["duty_mw"] = 0.0
        out["conduction_mw"] = 0.0
        out["radiation_mw"] = 0.0
        out["film_mw"] = 0.0
        return out
    m = re.search(
        r"wall duty: (-?[0-9.]+) MW leaving the gas into a (-?[0-9.]+) K wall",
        text,
    )
    if m:
        out["duty_mw"] = float(m.group(1))
        out["wall_k"] = float(m.group(2))
        out["has_duty_line"] = True
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
    out.setdefault("duty_source", "no duty line in the run log")
    return out


def analyse_leg(tag, runs_root, cases_dir):
    """everything measured about one leg."""
    run_dir = os.path.join(runs_root, tag)
    out = {"tag": tag, "ran": False}

    case_path = os.path.join(cases_dir, f"{tag}.toml")
    if not os.path.exists(case_path):
        out["why"] = "the case file is missing"
        return out

    out.update(read_wall_duty(os.path.join(run_dir, "run-log.txt")))

    fields = sorted(glob.glob(os.path.join(run_dir, "vtk-*.vtk")))
    if not fields:
        out["why"] = "no saved fields"
        return out

    # the flow falls into a strictly periodic orbit, so the measurement
    # is taken over WHOLE periods. a window covering part of a period
    # carries the wave's phase with it and is not comparable with
    # another leg's window.
    period = first = repeat = None
    if GATE is not None:
        period, first, repeat = GATE.detect_period(fields)
    if period is None and GATE is not None:
        # the orbit has not closed exactly, so the period is taken from
        # the observable's own recurrence, which the gate performs with
        # the same criterion, so the two instruments agree.
        fallback, strength = GATE.detect_period_by_signal(
            CF.read_case(case_path), fields
        )
        if fallback is not None:
            period = fallback
            out["period_confirmed_by_checksum"] = False
            out["period_signal_strength"] = strength
    if period is None:
        out["why"] = (
            "the saved fields show no repeat and no clear recurrence, so "
            "the flow's period cannot be identified"
        )
        out["fields_saved"] = len(fields)
        return out
    if "period_confirmed_by_checksum" not in out:
        out["period_confirmed_by_checksum"] = True
    out["period_fields"] = period
    out["fields_saved"] = len(fields)

    case = CF.read_case(case_path)
    groups = GATE.period_groups(fields, period)
    complete = [g for g in groups if len(g) == period]
    if not complete:
        out["why"] = "no whole period among the saved fields"
        return out
    out["whole_periods"] = len(complete)

    rows = []
    per_period = []
    for group in complete:
        measured = []
        for path in group:
            try:
                row = measure_field(case, path)
            except Exception as exc:  # pragma: no cover - a bad field is a real fault
                out.setdefault("field_errors", []).append(
                    "%s: %s" % (os.path.basename(path), exc)
                )
                continue
            if row:
                measured.append(row)
        if len(measured) == period:
            rows.extend(measured)
            per_period.append(measured)
    if not rows:
        out["why"] = "no field yielded a choked state"
        return out

    out["fields_measured"] = len(rows)
    out["field_rows"] = rows

    # the series that is characterised is the ONE PER PERIOD, because
    # that is the level the bracket compares.
    levels = []
    for measured in per_period:
        for key in (
            "choke_m_dot",
            "choke_p",
            "choke_t",
            "chamber_p0",
            "chamber_t0",
        ):
            values = [m[key] for m in measured if key in m]
            if values:
                levels.append((key, statistics.fmean(values)))

    labels = {
        "choke_m_dot": "choked_mass_flow_kg_s",
        "choke_p": "choke_pressure_pa",
        "choke_t": "choke_temperature_k",
        "chamber_p0": "chamber_total_pressure_pa",
        "chamber_t0": "chamber_total_temperature_k",
    }

    # the level of each period is what the bracket compares, so the
    # uncertainty is the spread of those levels.
    for key, label in labels.items():
        values = [level for name, level in levels if name == key]
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
        elif len(values) == 1 and measured:
            # only one whole period is available, so there is no spread
            # between periods. the period is split in half instead, and
            # the halves' disagreement bounds the bias the phase within
            # the period introduces. the same calibration the gate uses.
            first_half = [
                m[key] for m in measured[: len(measured) // 2] if key in m
            ]
            second_half = [
                m[key] for m in measured[len(measured) // 2 :] if key in m
            ]
            if first_half and second_half:
                half_mean = statistics.fmean(first_half)
                other_mean = statistics.fmean(second_half)
                disagreement = abs(half_mean - other_mean)
                stats["effective_sem"] = disagreement / 2.0 * math.sqrt(2.0)
                stats["effective_pct_sem"] = (
                    100.0 * stats["effective_sem"] / stats["mean"]
                    if stats["mean"]
                    else 0.0
                )
                stats["uncertainty_basis"] = "split half within one period"
                stats["half_disagreement"] = disagreement
        out[label] = stats
        out[f"{label}_behaviour"] = behaviour
        # the within-period variation is reported too, because it is the
        # oscillation the period averaging removes.
        within = []
        for measured in per_period:
            series = [m[key] for m in measured if key in m]
            if series:
                within.append(statistics.pstdev(series))
        if within:
            out[label]["within_period_sd"] = statistics.fmean(within)

    out["ran"] = True
    return out


def contrast(reference, other, key):
    """the difference between two legs, with its combined uncertainty."""
    a = reference.get(key)
    b = other.get(key)
    if not a or not b:
        return None
    ea = a.get("effective_sem", a.get("sem", 0.0))
    eb = b.get("effective_sem", b.get("sem", 0.0))
    delta = b["mean"] - a["mean"]
    combined = math.hypot(ea, eb)
    return {
        "delta": delta,
        "delta_pct": 100.0 * delta / a["mean"] if a["mean"] else 0.0,
        "combined_sem": combined,
        "sigma": delta / combined if combined > 0 else float("inf"),
    }


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
    cases_dir = args.cases_dir or locate(
        os.path.join(REPO_CASES, "engine-wall-bracket")
    )

    rows = [analyse_leg(tag, runs_root, cases_dir) for tag in LEGS]

    print(
        "%-18s %6s %13s %11s %13s %11s %10s"
        % ("leg", "fields", "mdot kg/s", "SEM %", "chamber T0 K", "SEM %", "duty MW")
    )
    for r in rows:
        if not r["ran"]:
            print("%-18s %6s  NOT RUN: %s" % (r["tag"], "-", r.get("why")))
            continue
        mf = r.get("choked_mass_flow_kg_s") or {}
        ct = r.get("chamber_total_temperature_k") or {}
        duty = r.get("duty_mw")
        print(
            "%-18s %6d %13.6f %11.4f %13.3f %11.4f %10s"
            % (
                r["tag"],
                r["fields_measured"],
                mf.get("mean", float("nan")),
                mf.get("effective_pct_sem", mf.get("pct_sem", 0.0)),
                ct.get("mean", float("nan")),
                ct.get("effective_pct_sem", ct.get("pct_sem", 0.0)),
                ("%.4f" % duty) if duty is not None else "?",
            )
        )
        beh = mf.get("_behaviour") or r.get("choked_mass_flow_kg_s_behaviour", {})
        if beh.get("adequate"):
            print(
                "    mass flow series: drift %+.4f %%  block-between %.5f vs within %.5f  lag-1 r %+.3f"
                % (
                    beh["drift_pct"],
                    beh["block_between_sd"],
                    beh["block_within_sd"],
                    beh["lag1_autocorrelation"],
                )
            )
        if r.get("has_duty_line"):
            print(
                "    wall: conduction %.4f MW, radiation %.4f MW, film %.4f MW"
                % (
                    r.get("conduction_mw", 0.0),
                    r.get("radiation_mw", 0.0),
                    r.get("film_mw", 0.0),
                )
            )
        else:
            print("    wall: %s" % r.get("duty_source"))

    by_tag = {r["tag"]: r for r in rows}
    print()
    print("contrasts against the adiabatic reference")
    reference = by_tag.get("wall-adiabatic")
    for tag in ("wall-radiating", "wall-filmcooled"):
        other = by_tag.get(tag)
        if not reference or not other or not reference["ran"] or not other["ran"]:
            continue
        for key, label in (
            ("choked_mass_flow_kg_s", "choked mass flow"),
            ("chamber_total_temperature_k", "chamber total temperature"),
            ("choke_pressure_pa", "choke pressure"),
        ):
            c = contrast(reference, other, key)
            if c:
                print(
                    "  %-18s %-26s %+0.6g  %+0.4f %%   combined SEM %.6g   %+.1f sigma"
                    % (
                        tag,
                        label,
                        c["delta"],
                        c["delta_pct"],
                        c["combined_sem"],
                        c["sigma"],
                    )
                )

    if args.show_fields:
        print()
        print("per-field series")
        for r in rows:
            if not r["ran"]:
                continue
            print(" ", r["tag"])
            for row in r.get("field_rows", []):
                print(
                    "    %-18s mdot %10.6f  p %12.1f  T0 %9.2f"
                    % (
                        row["field"],
                        row.get("choke_m_dot", float("nan")),
                        row.get("choke_p", float("nan")),
                        row.get("chamber_t0", float("nan")),
                    )
                )

    if args.json_path:
        payload = {"runs_root": runs_root, "cases": cases_dir, "legs": rows}
        with open(args.json_path, "w") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())

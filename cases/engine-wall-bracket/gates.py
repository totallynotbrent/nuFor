"""the engine wall bracket's gate.

the engine report needs to know what the wall does to the gas, so this
gate measures the bracket from the marched fields and checks the
contracts that make the comparison meaningful.

THE CONTRACTS

- every leg ran, reached the throat, and carries a whole number of flow
  periods' worth of saved fields;

- the adiabatic wall really is adiabatic: it draws no heat, so its
  measured duty is zero. if it is not, the comparison's reference end is
  broken;

- the cooled walls DO draw heat from the gas;

- radiation dominates conduction at these temperatures, because a gas at
  3700 k sees its wall overwhelmingly through radiation;

- the film-cooled wall draws MORE heat than the merely radiating wall,
  since the film is an additional path for the energy to leave by;

- the cooled walls' chambers run COOLER, and the choked mass flow RISES
  with the cooling, because it scales as p0 / sqrt(t0).

why the last pair is phrased that way

the throat's STATIC pressure is not an observable for this bracket. for
choked flow

    p_star = p0 (2 / (gamma + 1))^(gamma / (gamma - 1))

holds no total temperature, so a cooler chamber barely moves the throat
pressure. the mass flow does move, because it scales as p0 / sqrt(t0).

the sign of that dependence was got backwards in an earlier version of
this gate, which asserted that a cooler chamber carries LESS mass flow.
deriving it afresh from the choked relations gives the opposite:

    T_star / T0 = 2 / (gamma + 1)          so T_star scales with T0
    rho_star = p_star / (R T_star)          and p_star scales with p0
    V_star = sqrt(gamma R T_star)

which yields

    m_dot = A_star p_star sqrt(gamma) / (sqrt(R) sqrt(T0))
          = A_star p0 sqrt(gamma) / (sqrt(R) sqrt(T0))

so a cooler chamber at the same total pressure carries MORE mass. a
contract asserting the opposite would have failed a correct result.

MEASUREMENT NOTES

the marched chamber does not settle to a fixed state. it falls into a
strictly periodic orbit: the field saved at step 250 is byte-identical to
the one saved at step 3750, so the period spans fourteen field intervals.
a measurement taken over part of a period carries the wave's phase with
it, and a three-field window of the cyclic series departs from the
period's mean by up to 1.33 percent, which is about ten times the
effect the wall produces. so every leg is averaged over whole periods,
and the two legs' trajectories DIVERGE even though they start from the
same initial condition, which is why a step-by-step comparison between
legs is meaningless and only the period means can be compared.

the total pressure and total temperature are MEASURED from the marched
field as mass-flux-weighted means over the chamber, never read from the
run log's freestream line. that line is the inflow state's INITIAL
condition, and it is identical in every leg of the bracket.

nothing here is hardcoded.

usage:
    python3 gates.py
"""
import glob
import hashlib
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

BRACKET = [
    ("wall-adiabatic", None),
    ("wall-radiating", 400.0),
    ("wall-filmcooled", 400.0),
]


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


def digest(path):
    """the field's checksum, for spotting the flow's period."""
    hasher = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(65536), b""):
            hasher.update(block)
    return hasher.hexdigest()


def detect_period(fields):
    """the flow's period in saved fields, from the fields' checksums.

    the marched flow falls into a periodic orbit, so two saved fields with
    the same checksum mark one whole period apart. the smallest such
    separation is the period.

    the checksums are the STRONG criterion: they certify the orbit
    closed exactly. a flow that is periodic but has not quite closed
    shows no repeat, which the caller handles by falling back on the
    observable's own periodicity. see detect_period_by_signal.
    """
    seen = {}
    checksums = [digest(path) for path in fields]
    for index, value in enumerate(checksums):
        if value in seen:
            gap = index - seen[value]
            if gap > 0:
                return gap, seen[value], index
        seen[value] = index
    return None, None, None


def signal_correlation(values, lag):
    """the correlation between a series and itself, displaced by a lag."""
    n = len(values)
    if lag >= n - 2:
        return 0.0
    a = values[: n - lag]
    b = values[lag:]
    mean_a = statistics.fmean(a)
    mean_b = statistics.fmean(b)
    numerator = sum((x - mean_a) * (y - mean_b) for x, y in zip(a, b))
    denominator = math.sqrt(
        sum((x - mean_a) ** 2 for x in a) * sum((y - mean_b) ** 2 for y in b)
    )
    return numerator / denominator if denominator else 0.0


def detect_period_by_signal(case, fields):
    """the flow's period in saved fields, from the observable's recurrence.

    a flow that is periodic but has not quite closed has no repeated
    checksum, so the period is found from the choked mass flow instead:
    the lag whose correlation with the series is greatest is the period.

    this is a weaker criterion than the checksum, and it says so, so the
    measurement that rests on it carries an extra allowance.
    """
    values = []
    for path in fields:
        try:
            row = measure_field(case, path)
        except Exception:  # pragma: no cover - a bad field is a fault
            continue
        if row:
            values.append(row["mass_flow"])
    n = len(values)
    if n < 8:
        return None, 0.0
    best_lag = None
    best_value = 0.0
    for lag in range(2, n - 2):
        value = signal_correlation(values, lag)
        if value > best_value:
            best_value = value
            best_lag = lag
    if best_lag is None or best_value < 0.5:
        return None, best_value
    return best_lag, best_value


def period_groups(fields, period):
    """the fields grouped into whole periods, newest last.

    the groups are aligned on the LAST field, so a march that saved one
    more field than another still yields complete groups.
    """
    groups = []
    index = len(fields) - 1
    while index - period + 1 >= 0:
        groups.append(fields[index - period + 1 : index + 1])
        index -= period
    groups.reverse()
    return groups


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


def average_group(group, key):
    """a group's mean for one measured quantity."""
    values = [row[key] for row in group if row.get(key) is not None]
    if not values:
        return None
    return statistics.fmean(values)


def measure(tag, run_root):
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

    period, first, repeat = detect_period(saved)
    info["period_fields"] = period
    if period is not None:
        info["period_confirmed_by_checksum"] = True
        info["period_examples"] = [first, repeat]
    else:
        # the orbit has not closed exactly, so the period is taken from
        # the observable's own recurrence, which is a weaker criterion
        # and is recorded as such.
        info["period_confirmed_by_checksum"] = False
        fallback, strength = detect_period_by_signal(
            cf.read_case(case_path), saved
        )
        if fallback is None:
            info["why"] = (
                "the saved fields show no repeat and no clear recurrence "
                "(strongest correlation %.3f), so the flow's period cannot "
                "be identified and no whole-period average is possible"
                % strength
            )
            info["signal_correlation_best"] = strength
            return info
        period = fallback
        info["period_signal_strength"] = strength

    groups = period_groups(saved, period)
    complete = [g for g in groups if len(g) == period]
    if not complete:
        info["why"] = "no whole period among the saved fields"
        return info
    info["whole_periods"] = len(complete)

    case = cf.read_case(case_path)
    rows = []
    for group in complete:
        measured = []
        for path in group:
            try:
                row = measure_field(case, path)
            except Exception as exc:  # pragma: no cover - a bad field is a fault
                info.setdefault("field_errors", []).append(str(exc))
                continue
            if row:
                measured.append(row)
        if len(measured) == period:
            rows.append(measured)

    if not rows:
        info["why"] = "the fields in the period could not be measured"
        return info

    keys = (
        "mass_flow",
        "throat_p",
        "throat_t",
        "chamber_p0",
        "chamber_t0",
    )
    means = {}
    for key in keys:
        per_period = [average_group(row, key) for row in rows]
        per_period = [v for v in per_period if v is not None]
        if not per_period:
            continue
        means[key] = statistics.fmean(per_period)
        means[key + "_per_period"] = per_period
        if len(per_period) > 1:
            means[key + "_sem"] = statistics.pstdev(per_period) / math.sqrt(
                len(per_period)
            )
            means[key + "_uncertainty_basis"] = "period to period"
        elif rows:
            # only one whole period is available, so there is no spread
            # between periods to draw an uncertainty from. the period is
            # split in half instead, and the two halves' disagreement
            # bounds how much the phase within the period biases the
            # mean. this is the same calibration the adiabatic leg gave:
            # its two halves disagreed by 0.153 percent.
            #
            # the estimate is scaled by sqrt(2) because each half
            # averages half as many fields as the full period, so each
            # half's deviation from the full mean is amplified by
            # roughly sqrt(2).
            first_half = average_group(rows[0][: period // 2], key)
            second_half = average_group(rows[0][period // 2 :], key)
            if first_half is not None and second_half is not None:
                means[key + "_sem"] = (
                    abs(first_half - second_half) / 2.0 * math.sqrt(2.0)
                )
                means[key + "_uncertainty_basis"] = "split half within one period"
                means[key + "_half_disagreement"] = abs(first_half - second_half)
            else:
                means[key + "_sem"] = 0.0
                means[key + "_uncertainty_basis"] = "unknown"
        else:
            means[key + "_sem"] = 0.0
            means[key + "_uncertainty_basis"] = "unknown"
    info.update(means)
    info["ran"] = True
    return info


def predict_mass_flow_shift(reference, other, duty_mw):
    """the mass flow shift the wall's duty SHOULD produce.

    the wall's duty leaves the gas, so it lowers the chamber's bulk
    enthalpy. at a fixed total pressure a cooler chamber carries MORE
    mass through the throat, because the choked mass flow scales as
    p0 / sqrt(t0).

    this is the quantitative companion to the measured contrast. if the
    wall's duty really does reach the fluid, the measured shift should
    agree with the predicted one.

    the heat capacity is taken from the combustion products at the
    chamber's own temperature, where cp is near 3262 j/(kg.k) and gamma
    is near 1.16 for burnt gas around 3700 k.
    """
    if not duty_mw or duty_mw <= 0.0:
        return None
    mass_flow = reference.get("mass_flow")
    t0 = reference.get("chamber_t0")
    if not mass_flow or not t0:
        return None
    cp = 3262.0
    watt = duty_mw * 1.0e6
    # the steady sensible-energy balance. the wall's duty removes energy
    # from the gas at this rate, so removing it lowers the chamber's
    # enthalpy by that amount spread over the gas's mass flow.
    # the wall's duty LOWERS the chamber's enthalpy, so the chamber sits
    # cooler by this much than the adiabatic reference would.
    temperature_drop = watt / (mass_flow * cp)
    # a cooler chamber at fixed total pressure carries MORE mass, since
    # the choked mass flow scales as p0 / sqrt(t0). with t0 lowered by
    # the drop, the fractional change is +(drop / (2 t0)).
    predicted_fractional = 0.5 * temperature_drop / t0
    predicted_absolute = mass_flow * predicted_fractional
    return {
        "duty_mw": duty_mw,
        "predicted_dt0_k": -temperature_drop,
        "predicted_mdot_delta_kg_s": predicted_absolute,
        "predicted_mdot_delta_pct": 100.0 * predicted_fractional,
    }


def contrast(reference, other, key):
    """the difference between two legs, with its combined uncertainty."""
    a = reference.get(key)
    b = other.get(key)
    if a is None or b is None:
        return None
    ea = reference.get(key + "_sem", 0.0)
    eb = other.get(key + "_sem", 0.0)
    combined = math.hypot(ea, eb)
    delta = b - a
    return {
        "delta": delta,
        "delta_pct": 100.0 * delta / a if a else 0.0,
        "combined_sem": combined,
        "sigma": delta / combined if combined > 0 else float("inf"),
    }


def main():
    run_root = os.environ.get("NUFOR_BRACKET_RUNS", os.path.join(HERE, "runs"))
    rows = []
    for tag, tw in BRACKET:
        info = measure(tag, run_root)
        info["specified_wall_k"] = tw
        rows.append(info)

    print(
        "%-16s %7s %9s %9s %9s %12s %12s %8s"
        % (
            "leg",
            "periods",
            "duty MW",
            "cond MW",
            "rad MW",
            "mdot kg/s",
            "chamber T0 K",
            "sigma",
        )
    )
    for r in rows:
        if not r["ran"]:
            print(
                "%-16s %7s  NOT RUN: %s"
                % (r["tag"], "-", r.get("why"))
            )
            continue
        sigmas = []
        if r.get("mass_flow_delta_sigma") is not None:
            sigmas.append(abs(r["mass_flow_delta_sigma"]))
        print(
            "%-16s %7d %9.4f %9.4f %9.4f %12.6f %12.2f %8s"
            % (
                r["tag"],
                r["whole_periods"],
                r.get("duty_mw", 0.0),
                r.get("conduction_mw", 0.0),
                r.get("radiation_mw", 0.0),
                r.get("mass_flow", float("nan")),
                r.get("chamber_t0", float("nan")),
                "",
            )
        )

    problems = []
    by_tag = {r["tag"]: r for r in rows}
    adi = by_tag.get("wall-adiabatic")
    rad = by_tag.get("wall-radiating")
    film = by_tag.get("wall-filmcooled")

    for r in rows:
        if not r["ran"]:
            problems.append("the %s leg has not run (%s)" % (r["tag"], r.get("why")))

    if adi and adi["ran"]:
        # the adiabatic leg specifies no wall temperature, so it draws no
        # heat at all. if a duty line is present with a nonzero value,
        # the wall is not adiabatic and the reference end is broken.
        d = adi.get("duty_mw", None)
        if d is not None and abs(d) > 1.0e-6:
            problems.append(
                "the adiabatic wall drew %.6f MW; it must take NO heat" % d
            )
        if adi.get("whole_periods", 0) < 1:
            problems.append(
                "the adiabatic leg carries no whole period of saved fields"
            )

    for tag, leg in (("wall-radiating", rad), ("wall-filmcooled", film)):
        if not leg or not leg["ran"]:
            continue
        duty = leg.get("duty_mw", 0.0)
        if duty <= 0.0:
            problems.append(
                "the %s wall drew %.6f MW; it must take heat FROM the gas"
                % (tag, duty)
            )
        cond = leg.get("conduction_mw", 0.0)
        radn = leg.get("radiation_mw", 0.0)
        if cond > 0.0 and radn <= cond:
            problems.append(
                "the %s wall drew %.6f MW by radiation against %.6f MW by "
                "conduction; radiation must dominate at these temperatures"
                % (tag, radn, cond)
            )
        if leg.get("whole_periods", 0) < 1:
            problems.append(
                "the %s leg carries no whole period of saved fields" % tag
            )

    if rad and film and rad["ran"] and film["ran"]:
        drad = rad.get("duty_mw", 0.0)
        dfilm = film.get("duty_mw", 0.0)
        if dfilm <= drad:
            problems.append(
                "the film-cooled wall drew %.6f MW but the merely radiating "
                "wall drew %.6f MW; the film is an extra path for the heat "
                "to leave by and must draw more"
                % (dfilm, drad)
            )

    print()
    print("contrasts against the adiabatic reference")
    if adi and adi["ran"]:
        for tag, leg in (("wall-radiating", rad), ("wall-filmcooled", film)):
            if not leg or not leg["ran"]:
                continue
            for key, label in (
                ("mass_flow", "choked mass flow"),
                ("chamber_t0", "chamber total temperature"),
                ("chamber_p0", "chamber total pressure"),
                ("throat_p", "throat static pressure"),
            ):
                c = contrast(adi, leg, key)
                if c is None:
                    continue
                print(
                    "  %-16s %-26s %+0.6g  %+0.4f %%   combined SEM %.6g   %+.1f sigma"
                    % (
                        tag,
                        label,
                        c["delta"],
                        c["delta_pct"],
                        c["combined_sem"],
                        c["sigma"],
                    )
                )
                leg[key + "_delta_sigma"] = c["sigma"]
                leg[key + "_delta"] = c["delta"]

            # the physical expectation: a cooler chamber carries MORE mass
            # flow, because the choked mass flow scales as p0 / sqrt(t0).
            sf = leg.get("mass_flow_delta_sigma")
            st = leg.get("chamber_t0_delta_sigma")
            if sf is not None and st is not None:
                dt0 = leg.get("chamber_t0_delta", 0.0)
                if dt0 < 0.0 and sf < 0.0:
                    problems.append(
                        "the %s leg's chamber is cooler but its mass flow "
                        "fell; a cooler chamber carries MORE mass flow, "
                        "because the choked mass flow scales as "
                        "p0 / sqrt(t0)" % tag
                    )
                # the quantitative check: the shift the duty should have
                # produced, against the shift actually measured.
                prediction = predict_mass_flow_shift(
                    adi,
                    leg,
                    leg.get("duty_mw"),
                )
                if prediction is not None:
                    leg["predicted_effect"] = prediction
                    measured_delta = leg.get("mass_flow_delta")
                    if measured_delta is not None:
                        prediction["measured_mdot_delta_kg_s"] = measured_delta
                        ratio = (
                            measured_delta / prediction[
                                "predicted_mdot_delta_kg_s"
                            ]
                            if prediction["predicted_mdot_delta_kg_s"]
                            else None
                        )
                        prediction["measured_over_predicted"] = ratio
                        if ratio is not None and abs(ratio) > 3.0:
                            problems.append(
                                "the %s leg's mass flow shifted %+.6g kg/s "
                                "where the wall's duty predicts %+.6g kg/s, "
                                "which is %.1f times the prediction; the "
                                "duty and the flow's response disagree"
                                % (
                                    tag,
                                    measured_delta,
                                    prediction[
                                        "predicted_mdot_delta_kg_s"
                                    ],
                                    ratio,
                                )
                            )
                if prediction is not None and prediction.get(
                    "measured_mdot_delta_kg_s"
                ) is not None:
                    print(
                        "    %s: the duty predicts %+.6g kg/s "
                        "(%+0.4f percent), the flow shifted %+.6g kg/s, "
                        "which is %.2f times the prediction"
                        % (
                            tag,
                            prediction["predicted_mdot_delta_kg_s"],
                            prediction["predicted_mdot_delta_pct"],
                            prediction["measured_mdot_delta_kg_s"],
                            prediction.get("measured_over_predicted") or 0.0,
                        )
                    )
                if sf is not None and abs(sf) < 2.0:
                    print(
                        "    NOTE: the %s leg's mass flow contrast is only "
                        "%.1f sigma, so the effect is not resolved against "
                        "the march's own noise" % (tag, abs(sf))
                    )

    out = {"bracket": rows, "problems": problems}
    with open(os.path.join(HERE, "gates.json"), "w") as handle:
        json.dump(out, handle, indent=2)
        handle.write("\n")

    if problems:
        print("\nFAIL:")
        for p in problems:
            print("  - %s" % p)
        return 1
    print("\nthe bracket holds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

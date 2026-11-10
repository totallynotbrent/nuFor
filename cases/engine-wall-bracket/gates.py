"""the bracket's gate: what the campaign shows, judged against its contracts.

the bracket marches the same chamber under five wall treatments and
compares what each wall does to the gas. this module reads the
marched campaign, or the campaign's committed record, hands the
readings to `bracket.py`, and reports the verdict.

two inputs, one judgement. in the field the readings come from the
marched fields and the run logs; in continuous integration they come
from the record committed under `record/`. the contracts live in
`bracket.py` and are applied once, so the two cannot disagree about
what the campaign proved.

usage:

    python3 gates.py            # judge the marched campaign
    python3 gates.py --record   # judge the committed record
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
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "chamber-throat"))

import chamber_fields as cf  # noqa: E402

import bracket  # noqa: E402
import period_survey  # noqa: E402
import record  # noqa: E402
import wall_series  # noqa: E402


def canonical(value, places=9):
    """the verdict with its numbers rounded to a fixed number of places.

    the artefact travels between machines, and the arithmetic behind it
    lands on slightly different last bits depending on the platform's
    maths library. rounding to nine decimal places puts every quantity
    on a value the judgement actually claims, so the committed artefact
    compares the same way everywhere while a real change in the record
    still shows up.
    """
    if isinstance(value, float):
        if value != value or value in (float("inf"), float("-inf")):
            return value
        rounded = round(value, places)
        return 0.0 if rounded == 0.0 else rounded
    if isinstance(value, dict):
        return {key: canonical(item, places) for key, item in value.items()}
    if isinstance(value, list):
        return [canonical(item, places) for item in value]
    return value


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


def domain_mass(case, fields, nx, ny, origin, spacing):
    """the gas the chamber holds, integrated over the fluid cells.

    the peak density can climb for a benign reason: the gas cools against
    the wall and contracts. that is indistinguishable from mass being
    booked wrongly unless the chamber's total is watched too. this is the
    total, weighted the way the chamber weights its own columns.
    """

    mask = cf.solid_mask(case)
    dy = spacing[1]
    rho = fields["density"]
    total = 0.0
    for j in range(ny):
        rm = origin[1] + j * dy
        rp = rm + dy
        area = math.pi * (rp * rp - rm * rm)
        for i in range(nx):
            if mask[j * nx + i]:
                continue
            total += rho[j * nx + i] * area
    return total * spacing[0]


def measure_field(case, path):
    """one saved field's choked state and chamber totals."""
    fields, nx, ny, origin, spacing = cf.read_vtk(path)
    sonic = cf.sonic_point(case, fields, nx, ny, origin, spacing)
    if sonic is None:
        return None
    p0, t0 = chamber_totals(
        case, fields, nx, ny, origin, spacing, max(1, sonic["col"] - 1)
    )
    density = fields["density"]
    peak = max(range(len(density)), key=density.__getitem__)
    return {
        "path": path,
        "mass_flow": sonic["m_dot"],
        "throat_p": sonic["p"],
        "throat_t": sonic["T"],
        "throat_x": sonic["x"],
        "throat_basis": sonic["basis"],
        "chamber_p0": p0,
        "chamber_t0": t0,
        "domain_mass": domain_mass(case, fields, nx, ny, origin, spacing),
        "rho_max": density[peak],
        "t_at_rho_max": fields["temperature"][peak],
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

def field_index(path):
    """the march step a saved field came from, taken from its name."""
    stem = os.path.basename(path)
    if "." in stem:
        stem = stem.rsplit(".", 1)[0]
    digits = "".join(char for char in reversed(stem) if char.isdigit())
    return int(digits[::-1]) if digits else 0


def measure_one_leg(tag, run_root):
    """one leg of the bracket, read from the marched campaign."""
    run_dir = os.path.join(run_root, tag)
    case_path = os.path.join(HERE, "%s.toml" % tag)
    if not os.path.exists(case_path):
        return {"tag": tag, "ran": False, "why": "the case file is missing"}

    info = {"tag": tag, "ran": False}
    info.update(read_run_log(os.path.join(run_dir, "run-log.txt")))

    saved = sorted(glob.glob(os.path.join(run_dir, "vtk-*.vtk")))
    if not saved:
        info["why"] = "no saved fields"
        return info
    info["fields_saved"] = len(saved)

    rows = wall_series.read_duty_series(os.path.join(run_dir, "wall-duty.csv"))

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
        row["checksum"] = period_survey.digest(path)
        row["index"] = field_index(path)
        observations.append(row)

    info["density_track"] = bracket.measure_density_track(observations)
    return bracket.assemble_leg(info, rows, observations)


def measure_one_leg_from_record(tag):
    """one leg of the bracket, read from the committed record."""
    loaded = record.load_leg(tag)
    rows = wall_series.read_duty_series(record.record_paths(tag)["duty"])
    info = loaded["info"]
    info["density_track"] = bracket.measure_density_track(loaded["observations"])
    return bracket.assemble_leg(info, rows, loaded["observations"])


def judge(rows):
    """the bracket's judgement on some rows of legs.

    the rows carry each leg's measurements: its settled levels, its
    wall duty, and the saved fields behind them. the same judgement is
    made whether the rows came from marched fields or from the
    campaign's record.

    returns the report payload, the problems that make the bracket
    fail, and the cautions that bound what the numbers mean.
    """

    problems = []
    cautions = []
    by_tag = {r["tag"]: r for r in rows}

    # the settlement judgement is made on the wall's EFFECT: each leg's
    # departure from the adiabatic reference. the subtraction cancels the
    # chamber's own drift, which every leg carries and no leg owns, and
    # every leg is judged by the same rule, the reference included.
    cross = bracket.cross_leg_settlement.assess(rows)
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

    adi = by_tag.get(bracket.cross_leg_settlement.REFERENCE_TAG)
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

    for tag in bracket.COOLED_LEGS:
        leg = by_tag.get(tag)
        if not leg or not leg.get("ran"):
            continue
        departure = leg.get("mass_departure") or {}
        if departure and not departure.get("within_tolerance"):
            problems.append(
                "the %s wall departed from the adiabatic reference by %.3f "
                "percent of the chamber's mass, over the %.3f percent "
                "tolerance; the wall is booking mass the reference did not"
                % (
                    tag,
                    departure["delta_pct"],
                    departure["tolerance_pct"],
                ))
        duty = leg.get("duty_mw", 0.0)
        if duty <= 0.0:
            problems.append("the %s wall drew %.6f MW; it must take heat FROM the gas"
                            % (tag, duty))
        if tag not in bracket.COOLANT_BACKED_LEGS:
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
        for tag in bracket.COOLED_LEGS:
            leg = by_tag.get(tag)
            if not leg or not leg.get("ran"):
                continue
            for key, label in (
                ("mass_flow", "choked mass flow"),
                ("chamber_t0", "chamber total temperature"),
                ("chamber_p0", "chamber total pressure"),
                ("near_wall_t", "near-wall gas temperature"),
            ):
                c = bracket.paired_contrast(adi, leg, key)
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

            # the wall's own contribution to the chamber's mass. the
            # march's inventory drift is shared by every leg and is
            # reported in each leg's conservation block; what this
            # measures is the part the wall alone is responsible for.
            departure = bracket.paired_mass_departure(adi, leg)
            if departure is not None:
                leg["mass_departure"] = departure
                print("  %-17s %-26s %+0.6g  %+0.4f %%  of the chamber's mass" % (
                    tag, "domain mass", departure["delta"], departure["delta_pct"]))

            sf = leg.get("mass_flow_paired_contrast")
            if sf is not None:
                prediction = bracket.predict_choked_shift(adi, leg)
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
                scale = bracket.predict_duty_temperature_drop(adi, leg.get("duty_mw"))
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
    mode = bracket.measure_shared_mode(rows, "mass_flow")
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
    #
    # the survey reads the same way whichever side the measurements came
    # from: the checksum of each saved field and the mass flow read from
    # it. the field side takes its checksums from the files themselves;
    # the record side takes them from the record, which carries the same
    # digest beside each field's readings.
    per_leg = {}
    for row in rows:
        observations = row.get("observations") or []
        per_leg[row["tag"]] = {
            "digests": [
                observation["checksum"]
                for observation in observations
                if observation.get("checksum")
            ],
            "values": [
                observation["mass_flow"]
                for observation in observations
                if observation.get("mass_flow") is not None
            ],
        }
    survey = period_survey.survey(per_leg)

    print("period survey (report only; the gate does not rely on a period)")
    for tag, finding in sorted(survey.items()):
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

    return out, problems, cautions


def field_is_present(run_root):
    """whether the marched campaign is on disk here."""
    for tag, _ in bracket.BRACKET:
        run_dir = os.path.join(run_root, tag)
        if glob.glob(os.path.join(run_dir, "vtk-*.vtk")):
            return True
    return False


def read_the_campaign(run_root, from_record):
    """each leg's measurements, from the record or from the marched field."""
    rows = []
    marched = False
    for tag, wall_k in bracket.BRACKET:
        if from_record:
            info = measure_one_leg_from_record(tag)
        else:
            info = measure_one_leg(tag, run_root)
        info["specified_wall_k"] = wall_k
        rows.append(info)
        if info.get("ran"):
            marched = True
    return rows, marched


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--record",
        action="store_true",
        help="judge only the committed record, never the marched field",
    )
    parser.add_argument(
        "--no-cross-check",
        action="store_true",
        help="skip the field-versus-record cross check when both are present",
    )
    args = parser.parse_args(argv)

    run_root = os.environ.get("NUFOR_BRACKET_RUNS", os.path.join(HERE, "runs"))

    record_rows, _ = read_the_campaign(run_root, from_record=True)
    print("the engine wall bracket, judged against the committed record")
    out, problems, cautions = judge(record_rows)

    crossed = False
    if not args.record and not args.no_cross_check:
        if field_is_present(run_root):
            field_rows, _ = read_the_campaign(run_root, from_record=False)
            field_out, field_problems, field_cautions = judge(field_rows)
            crossed = True
            if field_problems:
                problems = problems + [
                    "the marched field fails where the record passes: %s" % problem
                    for problem in field_problems
                ]
            if field_cautions != cautions:
                problems.append(
                    "the marched field raises different cautions than the "
                    "record, so the record no longer describes the field"
                )
            disagreeing = [
                leg["tag"]
                for leg, field_leg in zip(out["bracket"], field_out["bracket"])
                if leg != field_leg
            ]
            if disagreeing:
                problems.append(
                    "the record and the marched field disagree on %s"
                    % ", ".join(disagreeing)
                )
            if not problems:
                print(
                    "the marched field agrees with the record on every leg"
                )
        else:
            print(
                "no marched field here, so the record stands on its own"
            )

    if not args.record:
        with open(os.path.join(HERE, "gates.json"), "w") as handle:
            json.dump(canonical(out), handle, indent=2)
            handle.write("\n")

    if cautions:
        print("\nCAUTIONS:")
        for caution in cautions:
            print("  - %s" % caution)

    if problems:
        print("\nFAIL:")
        for problem in problems:
            print("  - %s" % problem)
        return 1
    print("\nthe bracket holds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

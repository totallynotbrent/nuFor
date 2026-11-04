"""read the bracket's dense wall-duty series.

the bracket's run writes one `wall-duty.csv` per leg, sampled every
fifty steps independently of the field-saving interval. the series is
the bracket's densest observable, and unlike the saved fields it
covers the WHOLE run rather than one sample per interval.

WHY THE SERIES AND NOT THE FIELDS

the saved fields carry about a dozen samples across the run. that is
enough to see that the flow settled, but not enough to characterise
how well it settled: a dozen points give a crude trend estimate and a
poor variance estimate. the duty series gives sixty points at fifty
steps, sampled at the wall, which is where the bracket's physics
actually acts.

THE COLUMNS

    case, steps, time, heat_mw, conduction_mw, radiation_mw, film_mw,
    surface_t_k, perimeter_m

`heat_mw` is the NET power the wall took from the gas, as the solver
applied it. the three mode columns are MAGNITUDES of the individual
paths, so they do NOT sum to the net: the net is the mode sum projected
onto the wall's outward normal, and it is bounded jointly by the
step's admissibility limit. a decomposition that pretended the modes
added up would report energy that was never exchanged, so the gate
reads the net and the modes separately and checks the MODE ORDERING
rather than attempting a balance that does not exist.

this distinction matters for the adiabatic leg: its net is exactly zero
because the wall takes no heat at all, and its mode columns are zero
too, since an unwritten wall exchanges nothing by any path.
"""

import csv
import os

COLUMNS = (
    "steps",
    "time",
    "heat_mw",
    "conduction_mw",
    "radiation_mw",
    "film_mw",
    "surface_t_k",
    "perimeter_m",
)


def read_duty_series(path):
    """the leg's wall-duty series, oldest first.

    returns a list of dictionaries, one per row, with the numeric
    columns converted. a malformed row is skipped rather than fatal:
    a truncated final row is an artefact of a run stopping mid-write,
    not a reason to refuse the whole series.
    """
    if not os.path.exists(path):
        return []
    rows = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            row = {}
            usable = True
            for column in COLUMNS:
                value = raw.get(column)
                if value is None:
                    usable = False
                    break
                try:
                    row[column] = float(value)
                except (TypeError, ValueError):
                    usable = False
                    break
            if usable:
                rows.append(row)
    rows.sort(key=lambda r: r["steps"])
    return rows


def duty_summary(rows):
    """what the leg's wall did, summarized over the whole series.

    the figures are time-averaged over the series' own extent, which is
    the simplest honest summary when the series has not yet been
    divided into a settled window and an unsettled one. the gate's
    stationarity analysis then says how much those averages can be
    trusted.
    """
    if not rows:
        return {}

    def mean_of(column):
        values = [r[column] for r in rows]
        return sum(values) / len(values)

    net = mean_of("heat_mw")
    conduction = mean_of("conduction_mw")
    radiation = mean_of("radiation_mw")
    film = mean_of("film_mw")

    # the sign tells which way the energy flowed. the modes are
    # magnitudes, so the net's sign is the only statement about
    # direction the series can make.
    if net > 1.0e-9:
        direction = "the wall takes heat from the gas"
    elif net < -1.0e-9:
        direction = "the wall gives heat to the gas"
    else:
        direction = "the wall exchanges no heat"

    return {
        "duty_samples": len(rows),
        "duty_first_step": rows[0]["steps"],
        "duty_last_step": rows[-1]["steps"],
        "duty_mw": net,
        "conduction_mw": conduction,
        "radiation_mw": radiation,
        "film_mw": film,
        "surface_k": mean_of("surface_t_k"),
        "surface_perimeter_m": mean_of("perimeter_m"),
        "duty_direction": direction,
        "duty_net_sign": 1 if net > 1.0e-9 else (-1 if net < -1.0e-9 else 0),
    }


def mode_ordering(summary):
    """which mode leads the wall's loss, and by how much.

    the modes are magnitudes, so this ranks the PATHS rather than
    balancing them. the ranking is what the gate's contracts speak to:
    radiation leads on a bare wall, the coolant's film leads on a wall
    with a coolant behind it, and the liner's resistance shows up in
    the wall's surface temperature rather than in this ordering.
    """
    conduction = summary.get("conduction_mw", 0.0)
    radiation = summary.get("radiation_mw", 0.0)
    film = summary.get("film_mw", 0.0)
    modes = {
        "conduction": conduction,
        "radiation": radiation,
        "film": film,
    }
    leader = max(modes, key=lambda name: modes[name])
    total = sum(mode for mode in modes.values() if mode > 0.0)
    return {
        "leading_mode": leader,
        "leading_mode_mw": modes[leader],
        "mode_total_mw": total,
        "modes": modes,
    }

"""publishes the bracket's measured numbers into the page that quotes them.

THE PAGE SHOULD NOT CARRY HAND-TYPED MEASUREMENTS

every number this page quotes about the bracket was measured by the
gate, and the gate writes its judgement to `gates.json`. when the two
are kept apart by hand, they drift: a remarch moves the numbers and the
prose keeps quoting the old ones, and nothing catches it because the
prose has no way to disagree with itself.

so the measured numbers are published. this script reads the gate's
artefact and writes the page's marked regions, and the page keeps only
the prose around them. a number that leaves the artefact leaves the
page.

USAGE

    python3 publish_bracket_results.py <page> [--gates <path>]

the regions it owns are marked in the page with

    <!-- BEGIN TABLE <name> -->
    ...
    <!-- END TABLE <name> -->
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BEGIN = "<!-- BEGIN TABLE {name} -->"
END = "<!-- END TABLE {name} -->"

# the legs in the order the page introduces them.
ORDER = [
    "wall-adiabatic",
    "wall-radiating",
    "wall-filmcooled",
    "wall-regenerative",
    "wall-conjugate",
]


def legs_of(payload):
    """the artefact's legs, indexed by tag."""

    return {leg["tag"]: leg for leg in payload["bracket"]}


def mw(value):
    """a power in megawatts, to four decimals."""

    return "%.4f" % value


def kelvin(value):
    """a temperature, to one decimal."""

    return "%.1f" % value


def whole_kelvin(value):
    """a temperature rounded to the kelvin, for a wall the case states."""

    return "%g" % round(value)


def percent(value, digits=4):
    """a fraction as a percentage, to the requested digits."""

    return ("%." + str(digits) + "f") % value


def has_region(page, name):
    """whether the page declares a region with this name."""

    return BEGIN.format(name=name) in page


def replace_region(page, name, body):
    """swap one marked region for a new body."""

    begin = BEGIN.format(name=name)
    end = END.format(name=name)
    pattern = re.compile(re.escape(begin) + r".*?" + re.escape(end), re.DOTALL)
    if not pattern.search(page):
        return page, False
    body_text = begin + chr(10) + body.rstrip() + chr(10) + end
    return pattern.sub(lambda _: body_text, page), True


def duty_block(payload):
    """the wall duty each treatment drew, and how it divided."""

    legs = legs_of(payload)
    rows = []
    for tag in ORDER[1:]:
        leg = legs[tag]
        rows.append(
            "%-16s duty %s MW, conduction %s MW, radiation %s MW, film %s MW"
            % (
                tag.replace("wall-", "") + ":",
                mw(leg["duty_mw"]),
                mw(leg["conduction_mw"]),
                mw(leg["radiation_mw"]),
                mw(leg["film_mw"]),
            )
        )
    return "\n".join(rows)


def report_line_block(payload):
    """the film-cooled leg's duty as the march reports it."""

    legs = legs_of(payload)
    leg = legs["wall-filmcooled"]
    return (
        "wall duty: %s MW leaving the gas into a %s K wall\n"
        "  conduction %s MW\n"
        "  radiation  %s MW\n"
        "  film       %s MW"
        % (
            mw(leg["duty_mw"]),
            whole_kelvin(leg["wall_k"]) if leg.get("wall_k") else "400",
            mw(leg["conduction_mw"]),
            mw(leg["radiation_mw"]),
            mw(leg["film_mw"]),
        )
    )


def cooling_block(payload):
    """the cooling the duty implies against the cooling the gas registers."""

    legs = legs_of(payload)
    rows = []
    for tag in ORDER[1:]:
        leg = legs[tag]
        implied = leg["duty_scale_effect"]["scale_dt0_k"]
        measured = leg["near_wall_t_paired_contrast"]["delta"]
        rows.append(
            "%-16s implied %.1f K, measured %.1f K, difference %.1f K"
            % (
                tag.replace("wall-", "") + ":",
                abs(implied),
                abs(measured),
                abs(measured - implied),
            )
        )
    return "\n".join(rows)


def mass_flow_block(payload):
    """the choked mass flow each wall costs, and what the prediction said."""

    legs = legs_of(payload)
    rows = []
    for tag in ORDER[1:]:
        leg = legs[tag]
        contrast = leg["mass_flow_paired_contrast"]
        predicted = leg["predicted_effect"]
        rows.append(
            "%-16s measured %.3f percent, predicted %.3f percent, "
            "prediction error %.3f percent"
            % (
                tag.replace("wall-", "") + ":",
                abs(contrast["delta_pct"]),
                abs(predicted["predicted_mdot_delta_pct"]),
                predicted["measured_error_pct"],
            )
        )
    return "\n".join(rows)


def density_block(payload):
    """how each leg's peak density behaved, and why it moved."""

    legs = legs_of(payload)
    rows = []
    for tag in ORDER:
        leg = legs[tag]
        track = leg.get("density_track") or {}
        contraction = leg.get("contraction") or {}
        conservation = leg.get("conservation") or {}
        departure = leg.get("mass_departure") or {}
        if not track:
            rows.append("%-16s no density track recorded" % (tag.replace("wall-", "") + ":"))
            continue
        rows.append(
            "%-16s peak density climbs %.2f percent over the march "
            "(%.4f then %.4f percent per thousand steps), "
            "cooling %.1f K under the peak, chamber mass drifts %+.2f percent, "
            "wall's departure %+.3f percent of the chamber's mass"
            % (
                tag.replace("wall-", "") + ":",
                contraction.get("rho_max_climb_pct", 0.0),
                track.get("early_rate_pct_per_1000_steps", 0.0),
                track.get("late_rate_pct_per_1000_steps", 0.0),
                contraction.get("t_at_rho_max_drop_k", 0.0),
                conservation.get("inventory_drift_pct", 0.0),
                departure.get("delta_pct", 0.0),
            )
        )
    return "\n".join(rows)


def surface_block(payload):
    """the temperature the gas met at each wall."""

    legs = legs_of(payload)
    rows = []
    for tag in ORDER[1:]:
        leg = legs[tag]
        settled = ""
        if leg.get("surface_k_settled") is not None:
            settled = ", settled %.1f K" % leg["surface_k_settled"]
        rows.append(
            "%-16s duty %s MW into a %s K wall, surface %s K%s"
            % (
                tag.replace("wall-", "") + ":",
                mw(leg["duty_mw"]),
                whole_kelvin(leg["wall_k"]) if leg.get("wall_k") else "400",
                kelvin(leg["surface_k"]),
                settled,
            )
        )
    return "\n".join(rows)


def regenerative_surface_block(payload):
    """the regenerative wall's duty and the temperature its gas met."""

    legs = legs_of(payload)
    leg = legs["wall-regenerative"]
    return (
        "wall-regenerative  duty %s MW into a %s K wall\n"
        "                     (conduction %s MW, radiation %s MW, film %s MW)\n"
        "                   wall surface %s K, averaged over the wall's %s m of perimeter"
        % (
            mw(leg["duty_mw"]),
            whole_kelvin(leg["wall_k"]) if leg.get("wall_k") else "400",
            mw(leg["conduction_mw"]),
            mw(leg["radiation_mw"]),
            mw(leg["film_mw"]),
            kelvin(leg["surface_k"]),
            "%.3f" % leg["surface_perimeter_m"],
        )
    )


def settlement_block(payload):
    """the paired settlement each leg reached."""

    legs = legs_of(payload)
    rows = []
    for tag in ORDER[1:]:
        leg = legs[tag]
        contrast = leg.get("near_wall_t_paired_contrast") or {}
        rows.append(
            "%-16s near-wall departure %+.3f K against the adiabatic reference"
            % (tag.replace("wall-", "") + ":", contrast.get("delta", 0.0))
        )
    criterion = payload.get("settlement_criterion") or {}
    if criterion.get("qualified"):
        rows.append("")
        rows.append("The bracketing settlement criterion is met.")
    return "\n".join(rows)


def mode_block(payload):
    """the chamber's own motion, and the wall's imprint against it."""
    mode = payload.get("shared_mode") or {}
    levels = mode.get("levels") or {}
    ring = max((level.get("amplitude_pct", 0.0) for level in levels.values()), default=0.0)
    rows = []
    rows.append(
        "the chamber's mode: period %.0f steps, phase coherence %.4f, "
        "explains %.1f percent of the throat variation, "
        "worst cross-leg prediction %.2f"
        % (
            mode.get("period_steps", 0.0),
            mode.get("phase_coherence", 0.0),
            100.0 * mode.get("adjusted_explained_fraction", 0.0),
            mode.get("weakest_cross_prediction", 0.0),
        )
    )
    rows.append(
        "the mode's own amplitude on the throat: %.3f percent at its widest"
        % ring
    )
    legs = legs_of(payload)
    imprints = []
    for tag in ORDER[1:]:
        contrast = (legs[tag].get("mass_flow_paired_contrast") or {})
        imprints.append(abs(contrast.get("delta_pct", 0.0)))
    if imprints:
        rows.append(
            "the wall's imprint on the throat: at most %.3f percent, "
            "against the mode's %.3f percent"
            % (max(imprints), ring)
        )
    band = mode.get("period_band") or {}
    if band.get("lower_steps") is not None:
        rows.append(
            "the period's own tolerance: %.0f to %.0f steps hold the fit "
            "within %.1f percent of its best, and the ring stays between "
            "%.3f and %.3f percent across them"
            % (
                band.get("lower_steps", 0.0),
                band.get("upper_steps", 0.0),
                100.0 * band.get("tolerance_fraction", 0.0),
                band.get("worst_ring_percent_min", 0.0),
                band.get("worst_ring_percent_max", 0.0),
            )
        )
    return "\n".join(rows)


def screen_block(payload):
    """the settlement screen's own thresholds."""
    criterion = payload.get("settlement_criterion") or {}
    rows = []
    rows.append(
        "screens run: %d, over a window of %d saved fields (%d to %d)"
        % (
            criterion.get("comparisons", 0),
            criterion.get("window_length", 0),
            criterion.get("window_start", 0),
            criterion.get("window_end", 0),
        )
    )
    rows.append(
        "slope tolerance: %.2f standard errors, bonferroni-corrected over the family"
        % criterion.get("slope_tolerance", 0.0)
    )
    unsettled = criterion.get("secondary_unsettled") or []
    for note in unsettled:
        rows.append("unsettled secondary screen: " + note)
    return "\n".join(rows)


REGIONS = [
    ("report-line", report_line_block),
    ("duties", duty_block),
    ("cooling", cooling_block),
    ("mass-flow", mass_flow_block),
    ("density", density_block),
    ("surfaces", surface_block),
    ("regenerative-surface", regenerative_surface_block),
    ("settlement", settlement_block),
    ("mode", mode_block),
    ("screen", screen_block),
]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("page", help="the markdown page to write into")
    parser.add_argument("--gates", default=os.path.join(HERE, "gates.json"))
    args = parser.parse_args(argv)

    payload = json.load(open(args.gates))
    page = open(args.page).read()
    filled = []
    for name, builder in REGIONS:
        if not has_region(page, name):
            continue
        page, done = replace_region(page, name, builder(payload))
        if done:
            filled.append(name)
    if not filled:
        raise SystemExit("the page declares no published table region")
    open(args.page, "w").write(page)
    print("published %s into %s" % (", ".join(filled), args.page))


if __name__ == "__main__":
    main()

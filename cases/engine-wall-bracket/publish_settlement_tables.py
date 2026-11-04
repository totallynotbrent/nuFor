"""publishes the settlement tables from their JSON artefacts.

the bracket page quotes several measured tables: the certificate's
performance under drift, and the price of letting the screen choose its
own window. both come out of JSON artefacts produced by the harnesses.
this script rewrites the marked regions of the page from those files, so
a number on the page cannot drift away from the measurement that backs
it.

the regions are marked in the page with html comments of the form:

    <!-- BEGIN TABLE <name> -->
    ...
    <!-- END TABLE <name> -->

anything outside a marked region is left alone.

usage:

    python3 publish_settlement_tables.py <page> \\
        --power=<power.json> [--search=<search.json>]
"""

import argparse
import json
import os
import re


BEGIN = "<!-- BEGIN TABLE {name} -->"
END = "<!-- END TABLE {name} -->"


def percent(value):
    """a rate as a rounded percentage, or a dash when it is undefined."""
    if value is None:
        return "dash"
    return "{:.2f}".format(value * 100.0)


def powers_table(payload):
    """the detection table, one block per drift design."""
    blocks = []
    groups = {}
    for arm in payload["arms"]:
        groups.setdefault(arm["design"], []).append(arm)
    labels = {
        "whole_record": "Drift spread across the whole record",
        "confined": "Drift confined to the declared horizon",
    }
    for design in ("whole_record", "confined"):
        arms = groups.get(design)
        if not arms:
            continue
        blocks.append("")
        blocks.append(labels[design] + ":")
        blocks.append("")
        blocks.append("- drift, percent: detection among baseline, survival overall")
        for arm in arms:
            blocks.append(
                "- {:.3g}: {} percent, {} percent".format(
                    arm["amplitude_percent"],
                    percent(arm["detected_among_baseline"]),
                    percent(arm["certificate_survives_drift_rate"]),
                )
            )
    return "\n".join(block.strip("\n") for block in blocks)


def power_detail(payload):
    """the counts behind the detection table."""
    lines = []
    for arm in payload["arms"]:
        lines.append(
            "{} {:g} percent: baseline certified {}, drifted arm "
            "certified {}, withdrawn {}".format(
                arm["design"].replace("_", " "),
                arm["amplitude_percent"],
                arm["baseline_certificates"],
                arm["drift_arm_certificates"],
                arm["certificate_withdrawn"],
            )
        )
    return "\n".join(lines)


def search_table(payload):
    """the price of the window search, one block per drift design."""
    blocks = []
    groups = {}
    for cell_row in payload["cells"]:
        groups.setdefault(cell_row["design"], []).append(cell_row)
    labels = {
        "whole_record": "Drift spread across the whole record",
        "confined": "Drift confined to the declared horizon",
    }
    for design in ("whole_record", "confined"):
        cells = groups.get(design)
        if not cells:
            continue
        blocks.append("")
        blocks.append(labels[design] + ":")
        blocks.append("")
        blocks.append(
            "- drift, percent: search cost at matched width, "
            "unrestricted scan survival, registered survival"
        )
        for cell_row in cells:
            cost = cell_row["search_cost_at_matched_width"]
            blocks.append(
                "- {:g}: {} percentage points, {} percent, "
                "{} percent".format(
                    cell_row["amplitude_percent"],
                    percent(cost["miss_rate_difference"]),
                    percent(
                        cell_row["arms"]["scan_full"][
                            "certificate_survives_drift_rate"
                        ]
                    ),
                    percent(
                        cell_row["arms"]["registered"][
                            "certificate_survives_drift_rate"
                        ]
                    ),
                )
            )
    return "\n".join(block.strip("\n") for block in blocks)


def search_yield(payload):
    """what each arm accepts of a quiet record."""
    first = payload["cells"][0]["arms"]
    return (
        "On a quiet record the unrestricted scan certifies {} percent of "
        "surrogates, the scan capped at the horizon's width certifies {} "
        "percent, and the registered horizon certifies {} percent.".format(
            percent(first["scan_full"]["yield_on_quiet_record"]),
            percent(first["scan_width_matched"]["yield_on_quiet_record"]),
            percent(first["registered"]["yield_on_quiet_record"]),
        )
    )


BUILDERS = {
    "powers": powers_table,
    "power-detail": power_detail,
    "search-table": search_table,
    "search-yield": search_yield,
}


def replace_region(page, name, body):
    """swap one marked region for a new body."""
    begin = BEGIN.format(name=name)
    end = END.format(name=name)
    pattern = re.compile(
        re.escape(begin) + r".*?" + re.escape(end), re.DOTALL
    )
    if not pattern.search(page):
        raise SystemExit("the page has no table region named {}".format(name))
    return pattern.sub(begin + "\n" + body.rstrip() + "\n" + end, page)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("page", help="the markdown page to rewrite")
    parser.add_argument("--power", help="the detection-power JSON")
    parser.add_argument("--search", help="the search-multiplicity JSON")
    args = parser.parse_args(argv)
    payloads = {}
    if args.power:
        with open(args.power) as handle:
            payloads["power"] = json.load(handle)
    if args.search:
        with open(args.search) as handle:
            payloads["search"] = json.load(handle)

    with open(args.page) as handle:
        page = handle.read()

    if "power" in payloads:
        page = replace_region(page, "powers", powers_table(payloads["power"]))
        page = replace_region(
            page, "power-detail", power_detail(payloads["power"])
        )
    if "search" in payloads:
        page = replace_region(
            page, "search-table", search_table(payloads["search"])
        )
        page = replace_region(
            page, "search-yield", search_yield(payloads["search"])
        )

    with open(args.page, "w") as handle:
        handle.write(page)
    print("published the settlement tables into", os.path.basename(args.page))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

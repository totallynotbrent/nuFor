"""the bracket's record, and the tools that read and write it.

the bracket's five legs march the same chamber, so the interesting part
of the campaign is not one run but the DIFFERENCES between the runs.
those differences are expensive to recreate: the campaign costs hours
of wall clock and a field that no machine can rebuild on demand.

so the campaign is recorded once, in the raw form the judgement needs,
and committed. the judgement then runs against the RECORD in two places:

- in the field, where the full campaign lives and the gate reads the
  marched fields directly;
- in continuous integration, where the record stands in for the field.

that arrangement buys two things. the gate's contracts are checked
wherever anyone can run them, rather than only on the machine that
happened to march the campaign. and the published numbers have a single
source, so a page cannot drift away from the campaign without the gate
noticing.

the record holds RAW MEASUREMENTS, not conclusions. it stores the
per-field readings, the wall-duty series and the station breakdown for
each leg. the verdict is computed from those readings by the same code
that computed it in the field, so a stale or edited record fails the
gate instead of asserting itself.
"""

import csv
import json
import os


HERE = os.path.dirname(os.path.abspath(__file__))
RECORD_DIR = os.path.join(HERE, "record")


# the observables the record carries for every saved field. they are the
# ones the settlement judgement departs from and the ones the paired
# contrasts compare.
FIELD_KEYS = (
    "mass_flow",
    "chamber_p0",
    "chamber_t0",
    "near_wall_t",
    "domain_mass",
    "rho_max",
    "t_at_rho_max",
)


# the wall-duty columns the record carries for every march step.
DUTY_KEYS = (
    "steps",
    "time",
    "heat_mw",
    "conduction_mw",
    "radiation_mw",
    "film_mw",
    "surface_t_k",
    "perimeter_m",
)


def record_paths(tag):
    """the record's files for one leg."""
    leg_dir = os.path.join(RECORD_DIR, tag)
    return {
        "fields": os.path.join(leg_dir, "fields.csv"),
        "duty": os.path.join(leg_dir, "duty.csv"),
        "stations": os.path.join(leg_dir, "stations.csv"),
        "steady": os.path.join(leg_dir, "steadiness.csv"),
    }


def digest_column(index):
    """the field's checksum, as the record's spelling of it."""
    return "sha256:%s" % index


def write_csv(path, header, rows):
    """one csv file, written deterministically."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="\n") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def read_csv(path):
    """one csv file, read as a list of dictionaries."""
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def field_rows(observations, keys=FIELD_KEYS):
    """the per-field rows for a leg, ready for the csv."""
    rows = []
    for observation in observations:
        row = [observation["index"], observation.get("checksum", "")]
        for key in keys:
            value = observation.get(key)
            row.append("" if value is None else repr(float(value)))
        rows.append(row)
    return rows


def duty_table(duty_rows, keys=DUTY_KEYS):
    """the wall-duty rows for a leg, ready for the csv."""
    rows = []
    for row in duty_rows:
        cells = []
        for key in keys:
            value = row.get(key)
            if value in (None, ""):
                cells.append("")
            else:
                try:
                    cells.append(repr(float(value)))
                except (TypeError, ValueError):
                    cells.append(str(value))
        rows.append(cells)
    return rows


def station_table(stations):
    """the station rows for a leg, ready for the csv."""
    return [
        [repr(float(x0)), repr(float(x1)), repr(float(mw)), repr(float(perimeter))]
        for x0, x1, mw, perimeter in stations
    ]


def write_leg(tag, observations, duty_rows, stations):
    """one leg's raw measurements into the record."""
    paths = record_paths(tag)
    write_csv(
        paths["fields"],
        ("index",) + FIELD_KEYS,
        field_rows(observations),
    )
    write_csv(paths["duty"], DUTY_KEYS, duty_table(duty_rows))
    write_csv(
        paths["stations"],
        ("x0", "x1", "duty_mw", "perimeter_m"),
        station_table(stations),
    )


def load_leg(tag, specified_wall_k=None):
    """one leg's measurements, shaped the way the judge expects them."""
    paths = record_paths(tag)
    for key, path in sorted(paths.items()):
        if not os.path.exists(path):
            raise SystemExit(
                "the record for the %s leg is missing its %s at %s; "
                "the gate cannot judge a leg it cannot read"
                % (tag, key, path)
            )
    observations = []
    for row in read_csv(paths["fields"]):
        observation = {
            "path": "field-%05d.vtk" % int(row["index"]),
            "index": int(row["index"]),
        }
        checksum = row.get("checksum")
        if checksum:
            observation["checksum"] = checksum
        for key in FIELD_KEYS:
            value = row.get(key)
            if value in (None, ""):
                continue
            observation[key] = float(value)
        observations.append(observation)

    stations = []
    stations_path = paths["stations"]
    if os.path.exists(stations_path):
        for row in read_csv(stations_path):
            stations.append(
                (
                    float(row["x0"]),
                    float(row["x1"]),
                    float(row["duty_mw"]),
                    float(row["perimeter_m"]),
                )
            )

    info = {"tag": tag, "ran": True}
    scalars_path = os.path.join(os.path.dirname(paths["fields"]), "scalars.json")
    if os.path.exists(scalars_path):
        with open(scalars_path) as handle:
            scalars = json.load(handle)
        for key, value in scalars.items():
            if isinstance(value, (int, float)):
                info[key] = value

    return {
        "info": info,
        "observations": observations,
        "duty_rows": read_csv(paths["duty"]),
        "stations": stations,
        "specified_wall_k": specified_wall_k,
    }


def dump_json(path, payload):
    """one json file, written deterministically."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")

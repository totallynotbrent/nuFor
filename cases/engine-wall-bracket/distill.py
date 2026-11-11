"""records the bracket's campaign, so the gate can run without the fields.

the bracket's five legs march the same chamber for fifteen thousand
steps each, and the fields they save are far too large to keep in the
repository. what the gate needs from that campaign is small: the
per-field readings it judges the settlement on, the wall-duty series it
summarizes, and the station breakdown the page quotes.

this script distils the marched campaign into that small form and writes
it under `record/`. the record holds MEASUREMENTS, not conclusions: the
gate computes the verdict from the record with the same code it uses on
the field, so an edited or stale record fails the gate instead of
asserting itself.

usage:

    python3 distill.py [--runs <dir>] [--record <dir>]
"""

import argparse
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import figures  # noqa: E402
import gates  # noqa: E402
import period_survey  # noqa: E402
import record  # noqa: E402
import wall_series  # noqa: E402


def file_sha256(path):
    """a file's digest, so the record names the campaign it came from."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def field_observations(run_dir, tag):
    """the per-field readings the gate measures from a marched field.

    the readings come from the gate's own measurer, because that is the
    code the judgement has to agree with. the record therefore carries
    exactly what the gate would have measured.
    """
    case_path = os.path.join(HERE, "%s.toml" % tag)
    case = gates.cf.read_case(case_path)

    saved = sorted(
        name for name in os.listdir(run_dir) if name.startswith("vtk-")
    )
    observations = []
    for name in saved:
        path = os.path.join(run_dir, name)
        try:
            row = gates.measure_field(case, path)
        except Exception as exc:  # pragma: no cover - a bad field is a fault
            print("  %s: %s could not be measured: %s" % (tag, name, exc))
            continue
        if not row:
            continue
        try:
            row["near_wall_t"] = gates.near_wall_median(case, path, gates.cf)
        except Exception as exc:  # pragma: no cover
            print("  %s: %s has no near-wall reading: %s" % (tag, name, exc))
        row["index"] = int(name.split("-")[1].split(".")[0])
        row["checksum"] = period_survey.digest(path)
        observations.append(row)
    return observations


STEADY_KEYS = ("step", "t", "residual", "rho_max", "t_max")


def steady_rows(run_dir):
    """the march's own samples, so the page's steadiness claim is gated."""
    path = os.path.join(run_dir, "progress.json")
    if not os.path.exists(path):
        return []
    with open(path) as handle:
        data = json.load(handle)
    samples = data.get("samples") or []
    if not samples and isinstance(data.get("sample"), dict):
        samples = [data["sample"]]
    rows = []
    for sample in samples:
        rows.append([
            int(sample.get("step", 0)),
            repr(float(sample.get("t", 0.0))),
            repr(float(sample.get("residual", 0.0))),
            repr(float(sample.get("rho_max", 0.0))),
            repr(float(sample.get("t_max", 0.0))),
        ])
    return rows


def distil_leg(tag, run_root):
    """one leg's raw measurements, read from the marched campaign."""
    run_dir = os.path.join(run_root, tag)
    duty_path = os.path.join(run_dir, "wall-duty.csv")
    if not os.path.exists(duty_path):
        raise SystemExit("the %s leg has no duty series at %s" % (tag, duty_path))

    observations = field_observations(run_dir, tag)
    if not observations:
        raise SystemExit("the %s leg has no measurable fields" % tag)

    log_path = os.path.join(run_dir, "run-log.txt")
    log_text = open(log_path).read() if os.path.exists(log_path) else ""

    duty_rows = wall_series.read_duty_series(duty_path)
    stations = figures.parse_stations(log_text)
    scalars = gates.read_run_log(log_path)

    scalars["fields_saved"] = len(observations)
    scalars["duty_source_sha256"] = file_sha256(duty_path)

    return {
        "observations": observations,
        "duty_rows": duty_rows,
        "stations": stations,
        "scalars": scalars,
        "steady_rows": steady_rows(run_dir),
    }


def write_leg(tag, leg):
    """one leg's measurements into the record."""
    paths = record.record_paths(tag)
    os.makedirs(os.path.dirname(paths["fields"]), exist_ok=True)

    record.write_csv(
        paths["fields"],
        ("index", "checksum") + record.FIELD_KEYS,
        record.field_rows(leg["observations"]),
    )
    record.write_csv(
        paths["duty"],
        record.DUTY_KEYS,
        record.duty_table(leg["duty_rows"]),
    )
    record.write_csv(
        paths["stations"],
        ("x0", "x1", "duty_mw", "perimeter_m"),
        record.station_table(leg["stations"]),
    )
    record.write_csv(
        paths["steady"],
        STEADY_KEYS,
        leg["steady_rows"],
    )
    record.dump_json(
        os.path.join(os.path.dirname(paths["fields"]), "scalars.json"),
        leg["scalars"],
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", default=os.path.join(HERE, "runs"))
    parser.add_argument("--record", default=record.RECORD_DIR)
    args = parser.parse_args(argv)

    if not os.path.isdir(args.runs):
        raise SystemExit("there is no campaign at %s" % args.runs)

    manifest = {"legs": {}}
    for tag, _ in gates.bracket.BRACKET:
        print("reading the %s leg" % tag)
        leg = distil_leg(tag, args.runs)
        write_leg(tag, leg)
        manifest["legs"][tag] = {
            "case_sha256": file_sha256(os.path.join(HERE, "%s.toml" % tag)),
            "duty_source_sha256": leg["scalars"]["duty_source_sha256"],
            "fields_saved": leg["scalars"]["fields_saved"],
            "stations": len(leg["stations"]),
            "duty_rows": len(leg["duty_rows"]),
            "steadiness_source_sha256": file_sha256(
                os.path.join(args.runs, tag, "progress.json")
            ),
            "steadiness_samples": len(leg["steady_rows"]),
        }

    manifest["campaign"] = "the five legs of the engine wall bracket"
    manifest["measurements"] = "raw readings, judged by gates.py against this record"
    record.dump_json(os.path.join(args.record, "manifest.json"), manifest)

    print("recorded %d legs under %s" % (len(manifest["legs"]), args.record))
    return 0


if __name__ == "__main__":
    sys.exit(main())

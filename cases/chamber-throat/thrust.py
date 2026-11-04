"""thrust coefficients for the chamber-throat case.

reads the settled axial profile and computes the vacuum-style thrust
statement the roadmap promised for rung 9: the momentum-plus-pressure
thrust at the exit plane, referenced to the chamber pressure and throat
area as the thrust coefficient, and to gravity as the specific impulse.
mass flow is quoted at three stations so the gate can check the march
conserves it.

usage: python3 thrust.py <run-dir> <case.toml>
writes: <run-dir>/thrust.json
"""
import csv
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tomllib  # noqa: E402

G0 = 9.80665


def load_profile(run_dir):
    """rows of the settled profile csv."""
    path = os.path.join(run_dir, "settled", "profile.csv")
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    """compute and write the thrust record."""
    run_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "settled"
    )
    with open(sys.argv[2], "rb") as handle:
        case = tomllib.load(handle)
    rows = load_profile(run_dir)

    float_row = [{k: float(v) for k, v in row.items()} for row in rows]
    throat = min(float_row, key=lambda r: r["area_m2"])
    exit_plane = float_row[-1]
    chamber_cols = [r for r in float_row if r["x_m"] <= 0.075]
    chamber_p = sum(r["p_pa"] * r["area_m2"] for r in chamber_cols) / sum(
        r["area_m2"] for r in chamber_cols
    )
    chamber_mdot = sum(r["m_dot_kg_s"] for r in chamber_cols) / len(chamber_cols)

    pe = exit_plane["p_pa"]
    ae = exit_plane["area_m2"]
    ve = exit_plane["velocity_m_s"]
    mdot = exit_plane["m_dot_kg_s"]
    pa = case["boundaries"]["outflow_state"]["p"]
    at = throat["area_m2"]

    thrust = mdot * ve + (pe - pa) * ae
    cf = thrust / (chamber_p * at)
    isp = thrust / (mdot * G0)

    record = {
        "tag": "chamber-throat-thrust",
        "chamber_p0_pa": chamber_p,
        "throat_area_m2": at,
        "throat_x_m": throat["x_m"],
        "exit_area_m2": ae,
        "exit_p_pa": pe,
        "ambient_p_pa": pa,
        "exit_velocity_m_s": ve,
        "mass_flow_exit_kg_s": mdot,
        "mass_flow_throat_kg_s": throat["m_dot_kg_s"],
        "mass_flow_chamber_kg_s": chamber_mdot,
        "thrust_n": thrust,
        "thrust_coefficient": cf,
        "specific_impulse_s": isp,
        "checks": {
            "exit_supersonic": exit_plane["flux_mach"] > 1.0,
            "mass_conserved_across_stations": abs(mdot - throat["m_dot_kg_s"])
            / mdot
            < 0.15
            and abs(throat["m_dot_kg_s"] - chamber_mdot) / chamber_mdot < 0.6,
            "positive_thrust": thrust > 0.0,
            "reasonable_cf": 0.5 < cf < 2.5,
        },
    }
    record["passed"] = all(record["checks"].values())
    out = os.path.join(run_dir, "thrust.json")
    with open(out, "w") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    print(json.dumps(record, indent=2))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())

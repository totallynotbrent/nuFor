"""distils the chamber march into the record the docs and the gates read.

the march writes a snapshot every sample interval. this walks them, finds
where the flow reaches the speed of sound in each one (the mixed mach is
the criterion, and the throat's pressure is a property of that point,
not of the mesh's narrowest column), and writes:

    settled/history.csv        the sampled sonic and chamber states
    settled/settled.json       the settled state, its band and the profile
    settled/thermo_inputs.csv  what the rust audit needs to walk the
                               marched gas's own isentrope

usage:
    python3 distill.py <run dir>
"""
import argparse
import csv
import json
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import chamber_fields as cf  # noqa: E402

SPECIES = ["H2", "O2", "H2O", "OH", "H", "O", "HO2", "H2O2", "N2"]

# cea's equilibrium throat pressure for the case's chamber state
CEA_THROAT_P = 4.05352e6
CEA_THROAT_T = 3501.98

# the samples that count as the settled phase, the last ones of the march
SETTLED_SAMPLES = 8


def chamber_columns(case):
    """the columns that hold the chamber, up to where the duct converges."""
    dx = (case["x1"] - case["x0"]) / case["nx"]
    last = int((0.075 - case["x0"]) / dx)
    return range(5, max(6, last))


def analyse(path, case):
    """one snapshot: its sonic point, its chamber state and its gas."""
    fields, nx, ny, origin, spacing = cf.read_vtk(path)
    sonic = cf.sonic_point(case, fields, nx, ny, origin, spacing)
    columns = chamber_columns(case)
    chamber = cf.range_state(case, fields, nx, ny, origin, spacing, columns)
    gas = cf.composition(case, fields, nx, ny, origin, spacing, columns, SPECIES)
    stations = {}
    for col in range(nx):
        if col % 10 == 0 or col >= nx - 6:
            mix = cf.flux_composition(case, fields, nx, ny, origin, spacing, [col], SPECIES)
            if mix is not None:
                stations[col] = mix
    profile = []
    for col in range(nx):
        state = cf.station_state(case, fields, nx, ny, origin, spacing, col)
        if state is None:
            continue
        x = case["x0"] + (col + 0.5) * spacing[0]
        profile.append(
            {
                "x_m": x,
                "area_m2": state["area"],
                "p_pa": state["p"],
                "T_k": state["T"],
                "mixed_mach": state["mixed_mach"],
                "flux_mach": state["flux_mach"],
                "velocity_m_s": state["velocity"],
                "m_dot_kg_s": state["m_dot"],
            }
        )
    return sonic, chamber, gas, profile, stations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", help="the directory holding the march")
    parser.add_argument(
        "--case",
        default=os.path.join(HERE, "case.toml"),
        help="the case the march used",
    )
    parser.add_argument(
        "--out",
        default=os.path.join(HERE, "settled"),
        help="the directory for the distilled record",
    )
    args = parser.parse_args()

    case = cf.read_case(args.case)
    snapshots = sorted(
        name for name in os.listdir(args.run_dir) if name.startswith("vtk-") and name.endswith(".vtk")
    )
    if not snapshots:
        raise SystemExit(f"no snapshots in {args.run_dir}")
    step_of = lambda name: int(name.split("-")[1].split(".")[0])  # noqa: E731

    samples = []
    for name in snapshots:
        sonic, chamber, gas, profile, stations = analyse(os.path.join(args.run_dir, name), case)
        if sonic is None or chamber is None:
            continue
        samples.append(
            {
                "step": step_of(name),
                "sonic": sonic,
                "chamber": chamber,
                "gas": gas,
                "stations": stations,
                "profile": profile,
            }
        )
    samples.sort(key=lambda item: item["step"])
    if not samples:
        raise SystemExit("no snapshot carries a sonic point")

    settled = samples[-1]
    tail = samples[-SETTLED_SAMPLES:]
    band = lambda values: (statistics.mean(values), statistics.pstdev(values))  # noqa: E731

    sonic_p = [item["sonic"]["p"] for item in tail]
    sonic_t = [item["sonic"]["T"] for item in tail]
    sonic_mdot = [item["sonic"]["m_dot"] for item in tail]
    chamber_p = [item["chamber"]["p"] for item in tail]
    chamber_t = [item["chamber"]["T"] for item in tail]
    chamber_mdot = [item["chamber"]["m_dot"] for item in tail]

    out = args.out
    os.makedirs(out, exist_ok=True)

    with open(os.path.join(out, "history.csv"), "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "step",
                "sonic_col",
                "sonic_x_m",
                "sonic_p_pa",
                "sonic_T_k",
                "sonic_mdot_kg_s",
                "chamber_p_pa",
                "chamber_T_k",
                "chamber_mdot_kg_s",
            ]
        )
        for item in samples:
            writer.writerow(
                [
                    item["step"],
                    item["sonic"]["col"],
                    f"{item['sonic']['x']:.9e}",
                    f"{item['sonic']['p']:.6f}",
                    f"{item['sonic']['T']:.4f}",
                    f"{item['sonic']['m_dot']:.6f}",
                    f"{item['chamber']['p']:.6f}",
                    f"{item['chamber']['T']:.4f}",
                    f"{item['chamber']['m_dot']:.6f}",
                ]
            )

    record = {
        "case": "chamber-throat",
        "source": "the march distilled by cases/chamber-throat/distill.py",
        "step": settled["step"],
        "settled_samples": len(tail),
        "grid": [case["nx"], case["ny"]],
        "sonic": {
            "col": settled["sonic"]["col"],
            "x_m": settled["sonic"]["x"],
            "p_pa_mean": band(sonic_p)[0],
            "p_pa_sigma": band(sonic_p)[1],
            "T_k_mean": band(sonic_t)[0],
            "T_k_sigma": band(sonic_t)[1],
            "m_dot_kg_s_mean": band(sonic_mdot)[0],
            "m_dot_kg_s_sigma": band(sonic_mdot)[1],
        },
        "chamber": {
            "p_pa_mean": band(chamber_p)[0],
            "p_pa_sigma": band(chamber_p)[1],
            "T_k_mean": band(chamber_t)[0],
            "T_k_sigma": band(chamber_t)[1],
            "m_dot_kg_s_mean": band(chamber_mdot)[0],
            "m_dot_kg_s_sigma": band(chamber_mdot)[1],
        },
        "composition_mass_fraction": settled["gas"],
        "station_composition_mass_fraction": {str(col): mix for col, mix in sorted(settled["stations"].items())},
        "cea": {"throat_p_pa": CEA_THROAT_P, "throat_T_k": CEA_THROAT_T},
        "profile": settled["profile"],
    }
    with open(os.path.join(out, "settled.json"), "w") as handle:
        json.dump(record, handle, indent=1, sort_keys=True)
        handle.write("\n")


    with open(os.path.join(out, "stations.csv"), "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["col", "x_m"] + [f"y_{name}" for name in SPECIES])
        for col, mix in sorted(settled["stations"].items()):
            x = case["x0"] + (col + 0.5) * (case["x1"] - case["x0"]) / case["nx"]
            writer.writerow([col, f"{x:.9e}"] + [f"{mix[name]:.9f}" for name in SPECIES])

    with open(os.path.join(out, "profile.csv"), "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            ["x_m", "area_m2", "p_pa", "T_k", "mixed_mach", "flux_mach", "velocity_m_s", "m_dot_kg_s"]
        )
        for row in settled["profile"]:
            writer.writerow([
                f"{row['x_m']:.9e}",
                f"{row['area_m2']:.9e}",
                f"{row['p_pa']:.6f}",
                f"{row['T_k']:.4f}",
                f"{row['mixed_mach']:.6f}",
                f"{row['flux_mach']:.6f}",
                f"{row['velocity_m_s']:.6f}",
                f"{row['m_dot_kg_s']:.6f}",
            ])

    with open(os.path.join(out, "thermo_inputs.csv"), "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["# the chamber state and the sonic state the march measured"])
        writer.writerow(["chamber_p_pa", f"{band(chamber_p)[0]:.6f}"])
        writer.writerow(["chamber_T_k", f"{band(chamber_t)[0]:.4f}"])
        writer.writerow(["chamber_mach", f"{settled['chamber']['mixed_mach']:.6f}"])
        writer.writerow(["measured_sonic_x_m", f"{settled['sonic']['x']:.9e}"])
        writer.writerow(["measured_sonic_p_pa", f"{band(sonic_p)[0]:.6f}"])
        writer.writerow(["measured_sonic_T_k", f"{band(sonic_t)[0]:.4f}"])
        writer.writerow(["cea_throat_p_pa", f"{CEA_THROAT_P:.6f}"])
        for name in SPECIES:
            writer.writerow([f"y_{name}", f"{settled['gas'][name]:.9f}"])
        for row in settled["profile"]:
            writer.writerow([
                f"profile_{row['x_m']:.9e}_{row['p_pa']:.6f}_{row['T_k']:.4f}_{row['mixed_mach']:.6f}_{row['m_dot_kg_s']:.6f}".replace(",", "_"),
                "",
            ])

    print(
        f"samples {len(samples)}, settled sonic p {band(sonic_p)[0]/1e6:.4f} MPa "
        f"+/- {band(sonic_p)[1]/1e6:.4f}, chamber p {band(chamber_p)[0]/1e6:.4f} MPa"
    )


if __name__ == "__main__":
    main()

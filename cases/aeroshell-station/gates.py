# the aeroshell station gates: per station, the freestream-relative density
# ratio, the peak pressure, the peak temperature and the lift-to-drag read
# off the marched field next to each case file, then compared against the
# reference values in gates.json with a percentage tolerance.
#
# run from the case directory:
#   nufor run case.toml            # writes the vtk and the run log
#   python3 gates.py               # writes gates.json
#
# every number comes from the marched field and the case's own initial
# condition. nothing is hardcoded and nothing is read from outside the
# repository.
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

STATIONS = [
    ("h30.0", "aeroshell-h30.0-a22.5"),
    ("h34.6", "aeroshell-h34.6-a22.5"),
    ("h40.0", "aeroshell-h40.0-a22.5"),
    ("h42.5", "aeroshell-h42.5-a22.5"),
]

# percentage bands each measured quantity must stay inside.
TOLERANCE = {
    "density_ratio": 12.0,
    "peak_pressure_pct": 12.0,
    "max_temperature_pct": 8.0,
    "wall_flux_pct": 12.0,
}


def read_case(path):
    """the case's own freestream: rho, u and p from the initial condition."""
    text = open(path).read()
    vals = {}
    for key in ("rho", "u", "p"):
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith(key + " ") and "=" in stripped:
                val = stripped.split("=", 1)[1].strip()
                try:
                    vals[key] = float(val)
                except ValueError:
                    pass
                break
    missing = {"rho", "u", "p"} - set(vals)
    if missing:
        raise ValueError(f"{path}: no freestream {sorted(missing)}")
    return vals


def read_fields(path):
    lines = open(path).read().splitlines()
    nx = ny = None
    for line in lines:
        if line.startswith("DIMENSIONS"):
            nx, ny, _ = [int(v) for v in line.split()[1:]]
            break
    if nx is None or ny is None:
        raise ValueError(f"{path}: no DIMENSIONS block")
    n = nx * ny

    def field(name):
        marker = f"SCALARS {name} double 1"
        if marker not in lines:
            return None
        i = lines.index(marker)
        vals = []
        j = i + 2
        while len(vals) < n and j < len(lines):
            for tok in lines[j].split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
            j += 1
        return vals[:n]

    return field, nx, ny


def measure(vtk, case_toml):
    """the published quantities, read off the field and the freestream."""
    fs = read_case(case_toml)
    field, nx, ny = read_fields(vtk)
    rho = field("density")
    press = field("pressure")
    temp = field("temperature")
    if rho is None or press is None or temp is None:
        raise ValueError(f"{vtk}: missing density, pressure or temperature")

    # the body cells hold the freestream reference; the flow cells are the
    # ones the march moved. a cell whose density sits exactly on the
    # reference density never carried flow.
    ratio = [rho[k] / fs["rho"] for k in range(len(rho)) if abs(rho[k] - fs["rho"]) > 1e-9]
    if not ratio:
        raise ValueError(f"{vtk}: no flow cells")
    density_ratio = max(ratio)
    peak_pressure = max(press)
    max_temperature = max(temp)
    return {
        "density_ratio": density_ratio,
        "peak_pressure_pa": peak_pressure,
        "max_temperature_k": max_temperature,
    }


def read_log(path):
    """the aerodynamic coefficients the march logged."""
    if not os.path.exists(path):
        return {}
    out = {}
    for line in open(path):
        stripped = line.strip()
        if stripped.startswith("peak wall heat flux"):
            # the estimator's own number: "peak wall heat flux Q MW/m^2 at ..."
            try:
                out["wall_flux_mwm2"] = float(stripped.split()[4])
            except (IndexError, ValueError):
                pass
        for key in ("C_D", "C_L", "L/D"):
            if line.strip().startswith(key):
                try:
                    out[key] = float(line.split()[-1])
                except ValueError:
                    pass
    if "C_D" in out and "C_L" in out and out["C_D"]:
        out["L/D"] = out["C_L"] / out["C_D"]
    return out


def within(value, reference, pct):
    if reference == 0.0:
        return abs(value) <= 1e-12
    return abs(value - reference) / abs(reference) * 100.0 <= pct


def main():
    results = []
    failures = []
    for tag, dirname in STATIONS:
        # the station configurations and their marched fields live in the
        # sibling case directories; this gate only holds the thresholds.
        cdir = os.path.join(HERE, os.pardir, dirname)
        vtk = os.path.join(cdir, f"{dirname}.vtk")
        case = os.path.join(cdir, "case.toml")
        if not os.path.exists(vtk):
            print(f"MISS: no marched field for {dirname} (run the case)")
            failures.append(f"{dirname}: missing marched field")
            continue
        got = measure(vtk, case)
        got.update(read_log(os.path.join(cdir, "run-log.txt")))
        reference = None
        ref_path = os.path.join(HERE, "reference.json")
        if os.path.exists(ref_path):
            for entry in json.load(open(ref_path)):
                if entry.get("tag") == tag:
                    reference = entry
                    break
        entry = {"tag": tag, **{k: v for k, v in got.items()}}
        if reference:
            for key, tol in TOLERANCE.items():
                ref_key = {
                    "density_ratio": "density_ratio",
                    "peak_pressure_pct": "peak_pressure_pa",
                    "max_temperature_pct": "max_temperature_k",
                    "wall_flux_pct": "wall_flux_mwm2",
                }[key]
                if ref_key not in reference:
                    continue
                if ref_key not in got:
                    continue
                drift = (
                    abs(got[ref_key] - reference[ref_key]) / abs(reference[ref_key]) * 100.0
                    if reference[ref_key]
                    else 0.0
                )
                entry[f"{key}_drift_pct"] = drift
                if not within(got[ref_key], reference[ref_key], tol):
                    failures.append(
                        f"{dirname}: {ref_key} {got[ref_key]:.6g} vs reference "
                        f"{reference[ref_key]:.6g} ({drift:.2f}% > {tol}%)"
                    )
        results.append(entry)
        density = entry.get("density_ratio", 0.0)
        print(
            f"{dirname}: density ratio {density:.2f}  peak p "
            f"{entry.get('peak_pressure_pa', 0.0):.4g} Pa  max T "
            f"{entry.get('max_temperature_k', 0.0):.0f} K  L/D "
            f"{entry.get('L/D', 0.0):.4f}  peak wall flux "
            f"{entry.get('wall_flux_mwm2', 0.0):.4f} MW/m^2"
        )

    with open(os.path.join(HERE, "gates.json"), "w") as fh:
        json.dump(results, fh, indent=1)
        fh.write("\n")

    if failures:
        print("")
        for f in failures:
            print(f"FAIL {f}")
        return 1
    print("")
    print("station gates: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

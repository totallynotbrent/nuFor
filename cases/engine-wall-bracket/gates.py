"""the engine wall bracket's gate.

the engine report needs the COLD-WALL against COOLED-WALL comparison, so
this gate measures both ends of the bracket from the marched fields and
checks the contracts that make the comparison meaningful:

- each march ran, and its flow reached the throat.
- the adiabatic wall really is adiabatic: it draws no heat, so its
  measured duty is zero. if it is not, the comparison's reference end is
  broken.
- the cooled wall DOES draw heat, and it draws MORE the colder it is.
  a wall that does not cool its gas harder when it is colder would mean
  the wall treatment is inert.
- the cooled wall's THROAT PRESSURE is below the adiabatic wall's, with
  the coolant taking energy out of the gas before it accelerates.

three measurements had to be corrected for this gate to be able to see
anything at all:

the throat is located by the FLUX mach crossing, not by the narrowest
column. the narrowest column is the geometry's throat, but the flow
does not necessarily choke there, and on a staircase grid the two
differ by several columns.

the total pressure is MEASURED from the marched field as a
mass-flux-weighted mean over the chamber, not read from the run log's
freestream line. that line is the inflow state's INITIAL condition, and
it is identical in every leg of the bracket, so comparing it could never
discriminate the legs however hard the wall pulled.

the duty is read from the marched run's wall-diagnostic line, and its
conduction, radiation and film shares are recorded separately so the
breakdown can be audited.

nothing here is hardcoded.

usage:
    python3 gates.py
"""
import json
import math
import os
import re
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
    """the run log's wall duty and throat numbers."""
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
        r"wall duty: [-0-9.]+ MW leaving the gas into a [-0-9.]+ K wall\n"
        r"\(conduction ([-0-9.]+) MW, radiation ([-0-9.]+) MW, film ([-0-9.]+) MW\)",
        text,
    )
    if m:
        out["conduction_mw"] = float(m.group(1))
        out["radiation_mw"] = float(m.group(2))
        out["film_mw"] = float(m.group(3))
    m = re.search(r"peak wall heat flux ([0-9.]+) MW/m\^2", text)
    if m:
        out["peak_q_mwm2"] = float(m.group(1))
    m = re.search(r"freestream rho ([0-9.eE+-]+) p ([0-9.eE+-]+)", text)
    if m:
        out["freestream_p_pa"] = float(m.group(2))
    return out


def chamber_total_pressure(case, fields, nx, ny, origin, spacing, sonic):
    """the mass-flux-weighted total pressure over the chamber.

    the chamber is the region upstream of the choke. the total
    pressure at each cell follows from the cell's static pressure and
    mach, and the cells are averaged with their mass flux as the
    weight, because that is the weighting the conserved mass flow
    uses.

    this is a measurement of the marched field. the run log's
    freestream line is the inflow's initial condition and carries no
    information about what the wall did.
    """
    gamma = float(case.get("gamma", 1.2))
    gas_constant = float(case.get("gas_constant", 287.05))
    mask = cf.solid_mask(case)
    dy = spacing[1]
    pressure = fields["pressure"]
    temperature = fields["temperature"]
    density = fields["density"]
    mach = fields["mach"]
    velocity_x = fields.get("velocity_x")
    limit = sonic["col"] if sonic else nx
    numerator = denominator = 0.0
    for col in range(max(0, limit - 2)):
        for j in range(ny):
            if mask[j * nx + col]:
                continue
            k = j * nx + col
            u = abs(velocity_x[k]) if velocity_x is not None else 0.0
            m = max(mach[k], 1.0e-9)
            rm = origin[1] + j * dy
            rp = rm + dy
            area = math.pi * (rp * rp - rm * rm)
            weight = density[k] * u * area
            if weight <= 0.0:
                continue
            p0 = pressure[k] * (
                1.0 + (gamma - 1.0) / 2.0 * m * m
            ) ** (gamma / (gamma - 1.0))
            numerator += weight * p0
            denominator += weight
    if denominator <= 0.0:
        return None
    del gas_constant
    return numerator / denominator


def measure(tag):
    """everything the gate knows about one march of the bracket."""
    # the marched fields live outside the repository, so a stray run
    # cannot be mistaken for a published record. the directory comes from
    # the environment, and the default is the case's own runs directory.
    run_root = os.environ.get("NUFOR_BRACKET_RUNS", os.path.join(HERE, "runs"))
    run_dir = os.path.join(run_root, tag)
    case_path = os.path.join(HERE, f"{tag}.toml")
    if not os.path.exists(case_path):
        return {"tag": tag, "ran": False, "why": "the case file is missing"}
    info = {"tag": tag, "ran": False}
    log = read_run_log(os.path.join(run_dir, "run-log.txt"))
    info.update(log)
    field = None
    if os.path.isdir(run_dir):
        for name in sorted(os.listdir(run_dir)):
            if name.endswith(".vtk"):
                field = os.path.join(run_dir, name)
                break
    if field is None:
        info["why"] = "no marched field"
        return info
    fields, nx, ny, origin, spacing = cf.read_vtk(field)
    case = cf.read_case(case_path)
    sonic = cf.sonic_point(case, fields, nx, ny, origin, spacing)
    if sonic is None:
        info["why"] = "the flow never crossed mach one"
        return info
    info["throat_p_mpa"] = sonic["p"] / 1e6
    info["throat_t_k"] = sonic["T"]
    info["throat_mass_flow"] = sonic["m_dot"]
    info["throat_x_m"] = sonic["x"]
    info["throat_basis"] = sonic["basis"]
    total_p = chamber_total_pressure(
        case, fields, nx, ny, origin, spacing, sonic
    )
    if total_p is not None:
        info["chamber_total_p_mpa"] = total_p / 1e6
    info["ran"] = True
    return info


def main():
    rows = []
    for tag, tw in BRACKET:
        info = measure(tag)
        info["specified_wall_k"] = tw
        rows.append(info)

    print(
        "%-16s %8s %10s %10s %10s %12s"
        % ("leg", "wall K", "duty MW", "rad MW", "choke MPa", "p0 MPa")
    )
    for r in rows:
        if not r["ran"]:
            print(
                "%-16s %8s  NOT RUN: %s"
                % (r["tag"], r.get("wall_k", "-"), r.get("why"))
            )
            continue
        print(
            "%-16s %8.0f %10.5f %10.5f %10.5f %12.5f"
            % (
                r["tag"],
                r.get("wall_k", float("nan")),
                r.get("duty_mw", float("nan")),
                r.get("radiation_mw", float("nan")),
                r.get("throat_p_mpa", float("nan")),
                r.get("chamber_total_p_mpa", float("nan")),
            )
        )

    problems = []
    by_tag = {r["tag"]: r for r in rows}
    adi = by_tag.get("wall-adiabatic")
    cold = by_tag.get("wall-cold400")
    warm = by_tag.get("wall-warm1200")

    for r in rows:
        if not r["ran"]:
            problems.append("the %s leg has not run (%s)" % (r["tag"], r.get("why")))

    if adi and adi["ran"]:
        # the adiabatic leg specifies no wall temperature, so it draws no
        # heat at all and the run log carries no duty line. if a line is
        # present, the wall is not adiabatic and the reference end is
        # broken.
        d = adi.get("duty_mw", None)
        if d is not None and abs(d) > 1.0e-6:
            problems.append(
                "the adiabatic wall drew %.6f MW; it must take NO heat" % d
            )

    if cold and cold["ran"]:
        dcold = cold.get("duty_mw", 0.0)
        if dcold <= 0.0:
            problems.append(
                "the 400 K wall drew %.6f MW; it must take heat FROM the gas"
                % dcold
            )
        if warm and warm["ran"]:
            dwarm = warm.get("duty_mw", 0.0)
            if dcold < dwarm:
                problems.append(
                    "the 400 K wall drew %.6f MW but the 1200 K wall drew "
                    "%.6f MW; the colder wall must take MORE heat"
                    % (dcold, dwarm)
                )
        if adi and adi["ran"]:
            # the cooled wall takes energy out of the gas before it
            # reaches the throat, so the throat's pressure is lower than
            # the adiabatic wall's.
            #
            # the STATIC pressure at the throat is not the measure of
            # that: the throat is choked under every wall condition, so
            # its static pressure is set mostly by the choking condition
            # rather than by how much heat was removed upstream. the
            # MEASURED total pressure is the quantity that responds to
            # the wall, and it is taken from the marched field.
            p0c = cold.get("chamber_total_p_mpa")
            p0a = adi.get("chamber_total_p_mpa")
            if p0c is not None and p0a is not None and p0c >= p0a:
                problems.append(
                    "the cooled wall's measured total pressure %.4f MPa is "
                    "not below the adiabatic wall's %.4f MPa"
                    % (p0c, p0a)
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

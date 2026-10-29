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
- the cooled wall's throat pressure is BELOW the adiabatic wall's, with
  the coolant taking energy out of the gas before it accelerates.

every number comes from the marched field and the run log beside it.
nothing is hardcoded.

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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "chamber-throat"))
import chamber_fields as cf  # noqa: E402

BRACKET = [
    ("wall-adiabatic", None),
    ("wall-cold400", 400.0),
    ("wall-warm1200", 1200.0),
]


def read_run_log(path):
    """the run log's wall duty and throat numbers."""
    if not os.path.exists(path):
        return {}
    text = open(path).read()
    out = {}
    m = re.search(r"wall duty: (-?[0-9.]+) MW leaving the gas into a (-?[0-9.]+) K wall", text)
    if m:
        out["duty_mw"] = float(m.group(1))
        out["wall_k"] = float(m.group(2))
    m = re.search(r"peak wall heat flux ([0-9.]+) MW/m\^2", text)
    if m:
        out["peak_q_mwm2"] = float(m.group(1))
    return out


def measure(tag):
    """everything the gate knows about one march of the bracket."""
    run_dir = os.path.join(HERE, "runs", tag)
    case = os.path.join(HERE, f"{tag}.toml")
    if not os.path.exists(case):
        return {"tag": tag, "ran": False, "why": "the case file is missing"}
    info = {"tag": tag, "ran": False}
    log = read_run_log(os.path.join(run_dir, "run-log.txt"))
    info.update(log)
    settled = os.path.join(run_dir, "settled", "settled.json")
    if os.path.exists(settled):
        state = json.load(open(settled))
        info["sonic_p_mpa"] = state["sonic"]["p_pa_mean"] / 1e6
        info["sonic_t_k"] = state["sonic"]["T_k_mean"]
        info["mass_flow"] = state.get("mass_flow_kg_s")
        info["ran"] = True
    else:
        info["why"] = "no settled record"
    return info


def main():
    rows = []
    for tag, tw in BRACKET:
        info = measure(tag)
        info["specified_wall_k"] = tw
        rows.append(info)

    print("%-16s %10s %10s %12s %12s" % ("leg", "wall K", "duty MW", "sonic MPa", "sonic K"))
    for r in rows:
        if not r["ran"]:
            print("%-16s %10s  NOT RUN: %s" % (r["tag"], r.get("wall_k", "-"), r.get("why")))
            continue
        print(
            "%-16s %10.0f %10.5f %12.5f %12.1f"
            % (
                r["tag"],
                r.get("wall_k", float("nan")),
                r.get("duty_mw", float("nan")),
                r["sonic_p_mpa"],
                r["sonic_t_k"],
            )
        )

    problems = []

    adi = next((r for r in rows if r["tag"] == "wall-adiabatic"), None)
    cold = next((r for r in rows if r["tag"] == "wall-cold400"), None)
    warm = next((r for r in rows if r["tag"] == "wall-warm1200"), None)

    if adi and adi["ran"]:
        duty = abs(adi.get("duty_mw", 0.0))
        # the adiabatic wall draws no heat. the band's own allowance is
        # the numerical residue of the one-sided stencil.
        if duty > 1.0e-3:
            problems.append(
                "the adiabatic wall drew %.6f MW; it must draw none" % duty
            )
    else:
        problems.append("the adiabatic leg has not run")

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
            if cold["sonic_p_mpa"] >= adi["sonic_p_mpa"]:
                problems.append(
                    "the cooled wall's throat pressure %.5f MPa is not below "
                    "the adiabatic wall's %.5f MPa"
                    % (cold["sonic_p_mpa"], adi["sonic_p_mpa"])
                )
    else:
        problems.append("the cold wall leg has not run")

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

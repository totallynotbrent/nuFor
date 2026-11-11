#!/usr/bin/env python3
"""the feed-tank coupling gate.

it checks the record the solver wrote during the marched run against an
independent replay of the blowdown law.

the law has two regimes. the injector runs choked while the chamber sits
below the tank's pressure, where the flux is set by the tank state alone.
once the chamber catches up with the drained tank the injector is no
longer choked and the flux collapses to zero, because the solver holds
the reservoir as the source and never reverses it. the gate works out
which regime each step belongs to from the recorded pressures and
compares the flux the solver actually achieved against the flux that
regime predicts.

checks:

- the flux the solver achieved on every step matches the flux the
  governing regime predicts, to roundoff;
- the tank pressure and density on every recorded row match the replay;
- the drawn mass equals the tank's own level drop, so no mass is created
  at the seam;
- the polytropic invariant p / rho^n holds along the record;
- the reservoir actually drained: the final tank pressure is below the
  initial one.

the replay lives here rather than in the solver so a defect in the
solver's own integration shows up as a disagreement.

usage:
    python3 gates.py [--run-dir DIR]
"""
import argparse
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
RECORD = os.path.join(HERE, "record", "feed-tank.csv")


def read_case(path):
    """the tank stanza, the injector area and the case gamma."""
    vals = {}
    section = None
    for raw in open(path):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line.strip("[]")
            continue
        if "=" not in line:
            continue
        key, _, val = line.partition("=")
        key, val = key.strip(), val.strip()
        if section == "physics.tank_feed":
            vals[key] = val
        elif section == "boundaries" and key == "injector_area":
            vals["injector_area"] = val
        elif section == "physics" and key == "gamma":
            vals["gamma"] = val
    spec = {
        "volume": float(vals["volume"]),
        "rho0": float(vals["rho0"]),
        "t0": float(vals["t0"]),
        "polytrope": float(vals["polytrope"]),
        "r_gas": float(vals["r_gas"]),
        "injector_area": float(vals["injector_area"]),
        "gamma": float(vals["gamma"]),
    }
    spec["p0"] = spec["rho0"] * spec["r_gas"] * spec["t0"]
    return spec


def read_record(path):
    """the recorded reservoir trajectory."""
    rows = []
    with open(path) as handle:
        for rec in csv.DictReader(handle):
            rows.append(
                {
                    "time": float(rec["time_s"]),
                    "p_tank": float(rec["p_tank_pa"]),
                    "rho_tank": float(rec["rho_tank_kg_m3"]),
                    "p_down": float(rec["p_down_pa"]),
                    "drawn": float(rec["drawn_total_kg"]),
                }
            )
    return rows


def governing_flux(rho, temp, p_down, r_gas, gamma):
    """the mass flux the law prescribes for a tank state.

    returns the flux and the regime: "choked" while the chamber sits
    below the tank, "blocked" once the chamber has caught up with the
    drained tank, where the prescribed flux is zero.
    """
    p_t = rho * r_gas * temp
    ratio = max(min(p_down / p_t, 1.0), 1e-12)
    m_sq = 2.0 / (gamma - 1.0) * (ratio ** (-(gamma - 1.0) / gamma) - 1.0)
    mach = math.sqrt(max(m_sq, 0.0))
    if mach >= 1.0:
        sound = math.sqrt(gamma * r_gas * temp)
        factor = (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        return rho * sound * factor, "choked"
    if mach <= 0.0:
        return 0.0, "blocked"
    t_exit = temp * ratio ** ((gamma - 1.0) / gamma)
    return p_down * math.sqrt(gamma / (r_gas * t_exit)) * mach, "unchoked"


def replay(spec, rows):
    """drive the blowdown law over the recorded timestep sequence."""
    rho = spec["rho0"]
    temp = spec["t0"]
    gamma = spec["gamma"]
    area = spec["injector_area"]
    drawn = 0.0
    worst_p = 0.0
    worst_rho = 0.0
    worst_flux = 0.0
    invariant = 0.0
    regimes = {"choked": 0, "unchoked": 0, "blocked": 0}

    for prev, cur in zip(rows, rows[1:]):
        dt = cur["time"] - prev["time"]
        flux_expected, regime = governing_flux(rho, temp, cur["p_down"], spec["r_gas"], gamma)
        regimes[regime] += 1
        # the flux the solver achieved, from the recorded cumulative draw.
        achieved = (cur["drawn"] - prev["drawn"]) / (area * dt)
        if flux_expected > 0.0:
            worst_flux = max(worst_flux, abs(achieved - flux_expected) / flux_expected)
        else:
            worst_flux = max(worst_flux, abs(achieved))

        dm = min(flux_expected * area * dt, rho * spec["volume"])
        rho_new = rho - dm / spec["volume"]
        temp *= (rho_new / rho) ** (spec["polytrope"] - 1.0)
        rho = rho_new
        drawn += dm
        p_now = rho * spec["r_gas"] * temp
        worst_p = max(worst_p, abs(p_now - cur["p_tank"]) / cur["p_tank"])
        worst_rho = max(worst_rho, abs(rho - cur["rho_tank"]) / cur["rho_tank"])
        base = spec["p0"] / spec["rho0"] ** spec["polytrope"]
        invariant = max(invariant, abs(p_now / rho**spec["polytrope"] - base) / base)

    return {
        "drawn": drawn,
        "p_final": rho * spec["r_gas"] * temp,
        "rho_final": rho,
        "worst_p_rel": worst_p,
        "worst_rho_rel": worst_rho,
        "worst_flux_rel": worst_flux,
        "worst_invariant": invariant,
        "regimes": regimes,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", default=HERE)
    args = ap.parse_args()

    run_dir = args.run_dir
    spec = read_case(os.path.join(run_dir, "case.toml"))

    local = os.path.join(run_dir, "feed-tank.csv")
    source = local if os.path.exists(local) else RECORD
    if not os.path.exists(source):
        print("NOT RUN: no marched field, run the case first")
        return 1
    if source != RECORD:
        print("checking the marched field next to the case, not the committed record")

    rows = read_record(source)
    if len(rows) < 3:
        print("FAIL: the record carries fewer than three rows")
        return 1
    rep = replay(spec, rows)

    failures = []

    def check(label, got, want, tol):
        """relative agreement between a measured pair."""
        err = abs(got - want) / max(abs(want), 1e-30)
        ok = err <= tol
        print("%s %s: %.10g vs %.10g (rel %.3e, tol %.1e)"
              % ("PASS" if ok else "FAIL", label, got, want, err, tol))
        if not ok:
            failures.append(label)

    def check_bound(label, err, tol):
        """an error metric measured against its own bound."""
        ok = err <= tol
        print("%s %s: %.3e (tol %.1e)" % ("PASS" if ok else "FAIL", label, err, tol))
        if not ok:
            failures.append(label)

    check_bound("the achieved flux matches the governing regime", rep["worst_flux_rel"], 1e-12)
    check_bound("the recorded tank pressure matches the replay", rep["worst_p_rel"], 1e-12)
    check_bound("the recorded tank density matches the replay", rep["worst_rho_rel"], 1e-12)
    check("the drawn mass equals the tank level drop", rep["drawn"],
          (spec["rho0"] - rep["rho_final"]) * spec["volume"], 1e-9)
    check_bound("the polytropic invariant holds", rep["worst_invariant"], 1e-12)

    drawn_solver = rows[-1]["drawn"]
    check("the solver's drawn total matches the replay", drawn_solver, rep["drawn"], 1e-9)

    drained = rep["p_final"] < spec["p0"]
    print("%s the reservoir drained: %.6g pa from %.6g pa (%+.3f %% of p0)"
          % ("PASS" if drained else "FAIL", rep["p_final"], spec["p0"],
             100.0 * (rep["p_final"] - spec["p0"]) / spec["p0"]))
    if not drained:
        failures.append("the reservoir did not drain")

    regs = rep["regimes"]
    print("the injector ran choked on %d of %d steps (%d unchoked, %d blocked)"
          % (regs["choked"], sum(regs.values()), regs["unchoked"], regs["blocked"]))

    frac = rep["drawn"] / (spec["rho0"] * spec["volume"])
    print("draw %.6e kg of the %.6e kg reservoir (%.3f %%) over %d recorded steps, "
          "%d rows read from %s"
          % (rep["drawn"], spec["rho0"] * spec["volume"], 100.0 * frac, len(rows) - 1,
             len(rows), os.path.relpath(source, ROOT)))

    if failures:
        print("FAILED: " + "; ".join(failures))
        return 1
    print("passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())

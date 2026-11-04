#!/usr/bin/env python3.12
"""scan Cantera's ignition delay for the stoich methane-air parcel.

the chemistry page states that the same mixture ignites late at 1200 K and
not at all at 1000 K within a couple of milliseconds. those numbers come
from Cantera's own reactor, so this script re-runs that comparison and
prints the delays it finds, making the page's claim checkable rather than
merely asserted.

cantera 3.2 loads under python3.12 on this host and not under the older
interpreter, hence the version suffix in the name and in the shebang.

usage:
    env -u PYTHONPATH /usr/bin/python3.12 tools/scan_ignition_delays.py
"""

import argparse
import pathlib
import sys

CANTERA_SITE = "/home/ebarleb/.local/lib/python3.12/site-packages"


def build_parser() -> argparse.ArgumentParser:
    """describe the scan's options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--horizon",
        type=float,
        default=1.0,
        help="longest horizon to watch, in seconds",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=400.0,
        help="rise over the initial temperature that counts as ignition",
    )
    parser.add_argument(
        "--temperatures",
        default="1000,1200,1500,2000",
        help="comma-separated initial temperatures in kelvin",
    )
    return parser


def main(argv: list[str]) -> int:
    """run the scan and print a table of delays."""
    args = build_parser().parse_args(argv)

    if CANTERA_SITE not in sys.path:
        sys.path.insert(0, CANTERA_SITE)
    try:
        import cantera as ct
    except ImportError as exc:
        print(f"cannottera could not be imported: {exc}", file=sys.stderr)
        return 1

    mechanism = pathlib.Path(__file__).resolve().parents[1]
    mechanism = mechanism / "plans" / "combustion_research" / "sources" / "gri30.yaml"
    if not mechanism.is_file():
        print(f"mechanism not found: {mechanism}", file=sys.stderr)
        return 1

    gas = ct.Solution(str(mechanism))
    gas.TPX = 2000.0, 101325.0, "CH4:1, O2:2, N2:7.52"

    print(f"{'T0 (K)':>8} {'delay (s)':>14}  horizon reached")
    for text in args.temperatures.split(","):
        temperature = float(text)
        # the parcel starts unburned, at the temperature asked for and one
        # atmosphere, with the same mixture every time.
        gas.TPX = temperature, 101325.0, "CH4:1, O2:2, N2:7.52"

        # the plain reactor is closed, adiabatic-walled and fixed-volume,
        # which is exactly the parcel's situation. the copy argument keeps
        # the reactor and the gas sharing one object, so the unburned
        # composition stays intact.
        reactor = ct.Reactor(gas, clone=False)
        network = ct.ReactorNet([reactor])
        initial = gas.T
        delay = None
        step = min(args.horizon / 200.0, 5.0e-6)
        while network.time < args.horizon:
            try:
                # the integrator chooses its own steps but is brought to
                # each sample boundary, so the delay is read at a known
                # resolution rather than at the integrator's convenience.
                network.advance(network.time + step)
            except Exception as exc:  # noqa: BLE001
                print(f"{temperature:8.0f} {'failed':>14}  {exc}")
                break
            if reactor.T > initial + args.threshold:
                delay = network.time
                break
        if delay is None:
            print(f"{temperature:8.0f} {'no ignition':>14}  yes ({network.time:.4f} s)")
        else:
            print(f"{temperature:8.0f} {delay:14.6e}  no")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))

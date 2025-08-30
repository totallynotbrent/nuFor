#!/usr/bin/env python3
"""read a nufor hdf5 snapshot with h5py and print it as csv to stdout.

usage: python3 tools/export_h5.py snap.h5
columns: x,rho,m,e,u,p   (u and p derived with gamma=1.4)
"""

import h5py
import sys


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: export_h5.py <file.h5>")
    path = sys.argv[1]
    gamma = 1.4
    with h5py.File(path, "r") as f:
        x = f["centers"][:]
        rho = f["rho"][:]
        m = f["m"][:]
        e = f["e"][:]
    u = m / rho
    et = e / rho
    p = (gamma - 1.0) * rho * (et - 0.5 * u * u)
    print("x,rho,m,e,u,p")
    for i in range(len(x)):
        print(f"{x[i]},{rho[i]},{m[i]},{e[i]},{u[i]},{p[i]}")


if __name__ == "__main__":
    main()
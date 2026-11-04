"""drag coefficient from a marched park5-viscous vtk around the cylinder.

integrates the pressure traction over the fluid/solid interface: every
fluid cell whose centre lies just outside the body and touches a solid
cell contributes -p * n_hat . x_hat over its wetted face area, with n_hat
the outward body normal taken from the face midpoint. shear is neglected
at this reynolds number's leading order (casseau's ~1.3 is predominantly
pressure drag at m20).

    cd = f / (0.5 rho_inf u_inf^2 * 2 r)

run from the case directory:

    python3 drag_report.py cylinder-m20-park5-visc-full.vtk
"""
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

FREESTREAM = {"rho": 1.363e-5, "u": 6047.0}
RADIUS = 1.0
CX = CY = 0.0


def read_vtk(path):
    text = open(path).read()
    dim = re.search(r"DIMENSIONS (\d+) (\d+) (\d+)", text)
    ori = re.search(r"ORIGIN ([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)", text)
    sp = re.search(r"SPACING ([-\d.eE+]+) ([-\d.eE+]+) ([-\d.eE+]+)", text)
    nx, ny = int(dim.group(1)), int(dim.group(2))
    x0, y0 = float(ori.group(1)), float(ori.group(2))
    dx, dy = float(sp.group(1)), float(sp.group(2))
    fields = {}
    for m in re.finditer(r"SCALARS (\w+) double 1\nLOOKUP_TABLE default\n", text):
        name = m.group(1)
        start = m.end()
        vals = [float(v) for v in text[start:].split()[: nx * ny]]
        fields[name] = vals
    return nx, ny, x0, y0, dx, dy, fields


def main(vtk):
    nx, ny, x0, y0, dx, dy, f = read_vtk(vtk)
    p = f["pressure"]
    xs = [x0 + (i + 0.5) * dx for i in range(nx)]
    ys = [y0 + (j + 0.5) * dy for j in range(ny)]

    def inside(x, y):
        return (x - CX) ** 2 + (y - CY) ** 2 < RADIUS * RADIUS

    f_axial = 0.0
    n_iface = 0
    for j in range(ny):
        for i in range(nx):
            x, y = xs[i], ys[j]
            if inside(x, y):
                continue
            neigh = [(i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1)]
            for (ii, jj) in neigh:
                if ii < 0 or jj < 0 or ii >= nx or jj >= ny:
                    continue
                xn, yn = xs[ii], ys[jj]
                if inside(xn, yn):
                    k = j * nx + i
                    mx, my = 0.5 * (x + xn), 0.5 * (y + yn)
                    rl = math.hypot(mx - CX, my - CY) or 1e-30
                    nnx = (mx - CX) / rl
                    da = dy if ii != i else dx
                    f_axial += -p[k] * nnx * da
                    n_iface += 1
                    break
    q = 0.5 * FREESTREAM["rho"] * FREESTREAM["u"] ** 2
    cd = f_axial / (q * 2.0 * RADIUS)
    print(f"vtk={vtk}")
    print(f"interface_faces={n_iface}")
    print(f"axial_force_per_depth_n={f_axial:.6e}")
    print(f"dynamic_pressure_pa={q:.4e}")
    print(f"cd={cd:.4f}")


if __name__ == "__main__":
    main(sys.argv[1])

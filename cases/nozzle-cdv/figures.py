# nozzle cold-flow figures on the shared house template: the
# centerline mach vs the 1d area-mach solution for the three pressure
# ratios, and the wall-pressure distribution.
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from figtemplate import figure, save, ORANGE, PURPLE, GREEN, RED, GRAY


def read_fields(path):
    lines = open(path).read().split("\n")
    nx = ny = None
    for l in lines:
        if l.startswith("DIMENSIONS"):
            nx, ny, _ = [int(v) for v in l.split()[1:]]
            break

    def field(name):
        i = lines.index(f"SCALARS {name} double 1")
        vals = []
        j = i + 2
        while len(vals) < nx * ny and j < len(lines):
            vals += [float(v) for v in lines[j].split()]
            j += 1
        return vals[: nx * ny]

    return field, nx, ny


def read_1d():
    # regenerate the analytic reference from the case's own solver
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, here)
    import analytic_1d

    out = {}
    for pe, tag in ((0.89, "089"), (0.75, "075"), (0.16, "016")):
        sol, _kind, _crit = analytic_1d.solve_1d(pe)
        xs = [r[0] for r in sol]
        ms = [r[1] for r in sol]
        ps = [r[2] for r in sol]
        out[tag] = (xs, ms, ps)
    return out


def vtk_path(name):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), name)


cases = [
    ("089", "nozzle-cdv-p89.vtk", GREEN, "pe/pt = 0.89 (subsonic)"),
    ("075", "nozzle-cdv-p75.vtk", ORANGE, "pe/pt = 0.75 (normal shock)"),
    ("016", "nozzle-cdv-p16.vtk", RED, "pe/pt = 0.16 (supersonic)"),
]

one_d = read_1d()

fig, axes = figure(nrows=3, figsize=(7, 9))
for ax, (tag, vtk, color, label) in zip(axes, cases):
    xa, ma, _ = one_d[tag]
    ax.plot(xa, ma, color=GRAY, lw=1.2, ls="--", label="1d isentropic + shock")
    field, nx, ny = read_fields(vtk_path(vtk))
    m = field("mach")
    xs = [(i + 0.5) * 10.0 / nx for i in range(nx)]
    ax.plot(xs, [m[i] for i in range(nx)], color=color, lw=1.6, label="nuFor")
    ax.axvline(5.0, color=GRAY, lw=0.8, ls=":")
    ax.set_title(label)
    ax.set_ylabel("centerline Mach")
    ax.set_ylim(0, 3.2)
    ax.legend(loc="upper left")
    ax.annotate("throat", (5.06, 0.06), fontsize=9)
axes[-1].set_xlabel("x")
fig.tight_layout()
save(fig, "nozzle-cdv-mach", out_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "docs", "numerics", "2d"))

fig, ax = figure()
for tag, vtk, color, label in cases:
    field, nx, ny = read_fields(vtk_path(vtk))
    p = field("pressure")
    rho = field("density")
    dx = 10.0 / nx
    xs, ys = [], []
    for i in range(nx):
        x = (i + 0.5) * dx
        for j in range(ny - 1, -1, -1):
            k = j * nx + i
            if 1e-9 < rho[k] < 0.9:
                xs.append(x)
                ys.append(p[k])
                break
    ax.plot(xs, ys, color=color, lw=1.4, label=label + ", nuFor")
    xa, _, pa = one_d[tag]
    ax.plot(xa, pa, color=GRAY, lw=1.0, ls="--")
ax.axvline(5.0, color=GRAY, lw=0.8, ls=":")
ax.set_xlabel("x")
ax.set_ylabel("wall p/pt")
ax.set_title("wall pressure vs the 1d solution (gray dashed)")
ax.legend()
fig.tight_layout()
save(fig, "nozzle-cdv-wallp", out_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "docs", "numerics", "2d"))

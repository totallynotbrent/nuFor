# nozzle cold-flow rung figures: centerline mach vs the 1d
# area-mach solution for the three pressure ratios, plus the
# wall-pressure distribution. outputs into docs/research/.
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK = "#e8e8e8"; MUTE = "#999999"; CANVAS = "#141414"
GRID = "#2b2b2b"; AMBER = "#e0a458"; GREEN = "#7fb069"; RED = "#c1666b"

plt.rcParams.update({
    "figure.facecolor": CANVAS, "axes.facecolor": CANVAS,
    "text.color": INK, "axes.edgecolor": GRID, "axes.labelcolor": INK,
    "xtick.color": MUTE, "ytick.color": MUTE, "font.size": 11,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False,
})

def read_fields(path):
    lines = open(path).read().split('\n')
    nx = ny = None
    for l in lines:
        if l.startswith('DIMENSIONS'):
            nx, ny, _ = [int(v) for v in l.split()[1:]]
            break
    def field(name):
        i = lines.index(f"SCALARS {name} double 1")
        vals = []; j = i + 2
        while len(vals) < nx*ny and j < len(lines):
            vals += [float(v) for v in lines[j].split()]; j += 1
        return vals[:nx*ny]
    return field, nx, ny

def read_csv(path):
    xs, ms, ps = [], [], []
    for line in open(path):
        if line.startswith('#') or line.startswith('x'):
            continue
        x, m, p = line.strip().split(',')
        xs.append(float(x)); ms.append(float(m)); ps.append(float(p))
    return xs, ms, ps

cases = [
    ("p89", "p89_deep.vtk", "analytic_089.csv", GREEN,
     "pe/pt = 0.89 (subsonic)"),
    ("p75", "p75_deep.vtk", "analytic_075.csv", AMBER,
     "pe/pt = 0.75 (normal shock)"),
    ("p16", "p16_deep.vtk", "analytic_016.csv", RED,
     "pe/pt = 0.16 (fully supersonic)"),
]

fig, axes = plt.subplots(3, 1, figsize=(7.5, 10.5), sharex=True)
for ax, (tag, vtk, csv, color, label) in zip(axes, cases):
    xa, ma, pa = read_csv(csv)
    ax.plot(xa, ma, color=MUTE, lw=1.1, ls="--", label="1d isentropic + shock")
    field, nx, ny = read_fields(vtk)
    m = field('mach')
    xs = [(i + 0.5) * 10.0 / nx for i in range(nx)]
    ax.plot(xs, [m[i] for i in range(nx)], color=color, lw=1.4,
            label=label + ", nuFor")
    ax.axvline(5.0, color=GRID, lw=0.8)
    ax.set_ylabel("centerline Mach")
    ax.set_ylim(0, 3.2)
    ax.set_title(label, color=INK, fontsize=11)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.annotate("throat", (5.06, 0.06), color=MUTE, fontsize=8)
axes[-1].set_xlabel("x")
fig.tight_layout()
fig.savefig("../../docs/research/nozzle-cdv-mach.png", dpi=170,
            facecolor=CANVAS)
print("wrote mach figure")

fig, ax = plt.subplots(figsize=(7.5, 4.2))
for tag, vtk, csv, color, label in cases:
    field, nx, ny = read_fields(vtk)
    p = field('pressure'); rho = field('density')
    dx = 10.0 / nx
    xs, ys = [], []
    for i in range(nx):
        x = (i + 0.5) * dx
        for j in range(ny - 1, -1, -1):
            k = j * nx + i
            if rho[k] > 1e-9 and rho[k] < 0.9:
                xs.append(x); ys.append(p[k])
                break
    ax.plot(xs, ys, color=color, lw=1.2, label=label + ", nuFor")
    xa, ma, pa = read_csv(csv)
    ax.plot(xa, pa, color=MUTE, lw=1.0, ls="--")
ax.axvline(5.0, color=GRID, lw=0.8)
ax.set_xlabel("x"); ax.set_ylabel("wall p / pt")
ax.set_title("wall pressure vs the 1d solution (gray dashed)", color=INK,
             fontsize=11)
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("../../docs/research/nozzle-cdv-wallp.png", dpi=170,
            facecolor=CANVAS)
print("wrote wallp figure")

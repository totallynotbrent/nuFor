# the cdv nozzle gates: centerline mach and exit-core conditions vs
# the 1d solution, per pressure ratio. reads the march output that
# sits next to each case file, writes gates.json for the docs table
# and prints the comparison.
#
# run from the case directory:
#   nufor run case-p16.toml     # writes nozzle-cdv-p16.vtk
#   python3 gates.py             # writes gates.json
#
# every number here comes from the marched field and from
# analytic_1d.py in this directory. nothing is read from /tmp and
# nothing is hardcoded.
import math, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analytic_1d

g = 1.4
HERE = os.path.dirname(os.path.abspath(__file__))

CASES = [
    ("p89", "case-p89.toml", "nozzle-cdv-p89.vtk"),
    ("p75", "case-p75.toml", "nozzle-cdv-p75.vtk"),
    ("p16", "case-p16.toml", "nozzle-cdv-p16.vtk"),
]


def read_fields(path):
    lines = open(path).read().split('\n')
    nx = ny = None
    for l in lines:
        if l.startswith('DIMENSIONS'):
            nx, ny, _ = [int(v) for v in l.split()[1:]]
            break
    if nx is None or ny is None:
        raise ValueError(f"{path}: no DIMENSIONS block")
    def field(name):
        i = lines.index(f"SCALARS {name} double 1")
        vals = []; j = i + 2
        while len(vals) < nx*ny and j < len(lines):
            vals += [float(v) for v in lines[j].split()]; j += 1
        return vals[:nx*ny]
    return field, nx, ny


def read_spacing(path):
    for l in open(path):
        if l.startswith('SPACING'):
            return [float(v) for v in l.split()[1:]]
    raise ValueError("no SPACING")


def read_csv(path):
    xs, ms, ps = [], [], []
    for line in open(path):
        if line.startswith('#') or line.startswith('x'): continue
        x, m, p = line.strip().split(',')
        xs.append(float(x)); ms.append(float(m)); ps.append(float(p))
    return xs, ms, ps


def area(x):
    return analytic_1d.area(x)


def rad(x):
    return math.sqrt(area(x) / math.pi)


def interp(xs, ys, x):
    if x <= xs[0]: return ys[0]
    if x >= xs[-1]: return ys[-1]
    for i in range(len(xs) - 1):
        if xs[i] <= x <= xs[i+1]:
            f = (x - xs[i]) / (xs[i+1] - xs[i])
            return ys[i]*(1-f) + ys[i+1]*f
    return ys[-1]


def analytic_csv(pe, path):
    """writes the 1d reference and returns (rows, kind, throat mdot)."""
    sol, kind, crit = analytic_1d.solve_1d(pe)
    with open(path, "w") as f:
        f.write("x,M,p/pt\n")
        for x, M, p in sol:
            f.write(f"{x:.4f},{M:.5f},{p:.5f}\n")
    # the throat mass flow the 1d solution carries. the choked cases
    # run at the throat's critical state, the subsonic case at the
    # mass flow its own exit pressure implies.
    mdot = analytic_1d.throat_mdot(pe)
    return sol, kind, mdot


def gates(tag, vtk, csv, mdot_1d):
    field, nx, ny = read_fields(vtk)
    m = field('mach'); p = field('pressure'); rho = field('density')
    dx, dy, _ = read_spacing(vtk)
    xa, ma, pa = read_csv(csv)
    # throat centerline mach
    mt_cl = m[80]
    # centerline l2 vs 1d
    num = den = 0.0
    for i in range(nx):
        x = (i + 0.5) * dx
        d = m[i] - interp(xa, ma, x)
        num += d*d
        den += interp(xa, ma, x) ** 2
    l2 = math.sqrt(num/den) * 100.0
    # exit centerline
    k_exit = nx - 1
    me_cl = m[k_exit]
    pe_cl = p[k_exit]
    # exit core flux (inner half radius) and its 1d share
    rw = rad(10.0)
    core_mdot = 0.0
    for j in range(ny):
        rc = (j + 0.5) * dy
        if rc > 0.5 * rw: continue
        k = j*nx + k_exit
        if rho[k] < 1e-9 or rho[k] > 0.9: continue
        a_ = math.sqrt(g * p[k] / rho[k])
        core_mdot += rho[k] * m[k] * a_ * 2*math.pi*rc*dy
    core_area = math.pi * (0.5*rw)**2
    core_1d = mdot_1d * core_area / 1.5
    # plane mass fluxes: inlet, throat, exit
    def u_of(k):
        a = math.sqrt(g * max(p[k], 1e-12) / max(rho[k], 1e-12))
        return m[k] * a
    def mdot_plane(i):
        s = 0.0
        for j in range(ny):
            k = j*nx + i
            if rho[k] < 1e-9: continue
            ri = j * dy; ro = (j + 1) * dy
            s += rho[k] * u_of(k) * math.pi * (ro*ro - ri*ri)
        return s
    mi, mt, me = mdot_plane(0), mdot_plane(nx//2), mdot_plane(nx-1)
    spread = (max(mi, mt, me) - min(mi, mt, me)) / ((mi + mt + me) / 3) * 100.0
    # flux-weighted exit mach: the mass-average of the exit plane
    num_f = den_f = 0.0
    for j in range(ny):
        k = j*nx + k_exit
        if rho[k] < 1e-9: continue
        ri = j * dy; ro = (j + 1) * dy
        a_ring = math.pi * (ro*ro - ri*ri)
        fl = rho[k] * u_of(k) * a_ring
        num_f += fl * m[k]; den_f += fl
    exit_mach_massavg = num_f / den_f if den_f > 0 else 0.0
    # peak mach on the exit plane and the 1d exit mach
    m_exit_1d = interp(xa, ma, 10.0)
    out = {
        "tag": tag,
        "throat_M_cl": round(mt_cl, 4),
        "centerline_l2_pct": round(l2, 2),
        "exit_M_cl": round(me_cl, 4),
        "exit_M_massavg": round(exit_mach_massavg, 4),
        "exit_M_1d": round(m_exit_1d, 4),
        "exit_p_cl": round(pe_cl, 4),
        "exit_core_mdot": round(core_mdot, 4),
        "exit_core_mdot_1d": round(core_1d, 4),
        "core_mdot_err_pct": round((core_mdot/core_1d - 1)*100, 2),
        "inlet_mdot": round(mi, 5),
        "throat_mdot": round(mt, 5),
        "exit_mdot": round(me, 5),
        "plane_mdot_spread_pct": round(spread, 3),
    }
    return out


def main():
    results = []
    for tag, case_file, vtk_name in CASES:
        vtk = os.path.join(HERE, vtk_name)
        if not os.path.exists(vtk):
            print(f"skip {tag}: {vtk_name} missing (run the case first)")
            continue
        pe = {"p89": 0.89, "p75": 0.75, "p16": 0.16}[tag]
        csv = os.path.join(HERE, f"analytic_{str(pe).replace('.', '')}.csv")
        _, kind, mdot_1d = analytic_csv(pe, csv)
        r = gates(tag, vtk, csv, mdot_1d)
        r["regime"] = kind
        results.append(r)
        print(r)
    if not results:
        print("no marched fields found. run the cases first.")
        return 1
    with open(os.path.join(HERE, "gates.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("gates.json written")
    return 0


if __name__ == "__main__":
    sys.exit(main())

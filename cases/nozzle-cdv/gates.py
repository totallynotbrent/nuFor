# the cdv nozzle gates: centerline mach and exit-core conditions vs
# the 1d solution, per pressure ratio. writes gates.json for the docs
# table and prints the comparison.
import math, json

g = 1.4

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
        if line.startswith('#') or line.startswith('x'): continue
        x, m, p = line.strip().split(',')
        xs.append(float(x)); ms.append(float(m)); ps.append(float(p))
    return xs, ms, ps

def area(x):
    if x < 5.0:
        return 1.75 - 0.75 * math.cos((0.2 * x - 1.0) * math.pi)
    return 1.25 - 0.25 * math.cos((0.2 * x - 1.0) * math.pi)

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

def gates(tag, vtk, csv):
    field, nx, ny = read_fields(vtk)
    m = field('mach'); p = field('pressure'); rho = field('density')
    dx, dy = 10.0/nx, None
    for l in open(vtk):
        if l.startswith('SPACING'):
            dx, dy, _ = [float(v) for v in l.split()[1:]]
            break
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
    mdot_1d = float(open(f"mdot1d_{tag}.txt").read())
    core_1d = mdot_1d * core_area / 1.5
    out = {
        "tag": tag,
        "throat_M_cl": round(mt_cl, 4),
        "centerline_l2_pct": round(l2, 2),
        "exit_M_cl": round(me_cl, 4),
        "exit_p_cl": round(pe_cl, 4),
        "exit_core_mdot": round(core_mdot, 4),
        "exit_core_mdot_1d": round(core_1d, 4),
        "core_mdot_err_pct": round((core_mdot/core_1d - 1)*100, 2),
    }
    return out

if __name__ == "__main__":
    import sys
    mdot1d = {"p89": 0.3815, "p75": 0.3953, "p16": 0.3953}
    results = []
    for tag in ["p89", "p75", "p16"]:
        csv = f"analytic_0{tag[1:]}.csv"
        open(f"mdot1d_{tag}.txt", "w").write(f"{mdot1d[tag]}")
        r = gates(tag, f"/tmp/noz_runs/{tag}_deep.vtk", csv)
        results.append(r)
        print(r)
    json.dump(results, open("gates.json", "w"), indent=1)
    print("gates.json written")

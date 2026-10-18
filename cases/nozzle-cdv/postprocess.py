# cdv nozzle gates: centerline mach vs the 1d area-mach solution, wall
# pressure distribution, mass-flux conservation, exit conditions.
import math, sys

g = 1.4
def area(x):
    if x < 5.0:
        return 1.75 - 0.75 * math.cos((0.2 * x - 1.0) * math.pi)
    return 1.25 - 0.25 * math.cos((0.2 * x - 1.0) * math.pi)

def mach_from_area_ratio(ar, supersonic):
    # inverse area-mach: ar = A/A* (ar >= 1)
    if ar < 1.0 - 1e-12:
        return float("nan")
    f = lambda M: (1.0/M)*((2/(g+1))*(1+(g-1)/2*M*M))**((g+1)/(2*(g-1))) - ar
    lo, hi = (1.0, 12.0) if supersonic else (1e-6, 1.0)
    for _ in range(100):
        mid = 0.5*(lo+hi)
        if (f(mid) < 0) == (lo < 1.0):  # f increasing on subsonic branch? no:
            # subsonic branch: f decreasing from +inf at M->0 to 0 at M=1.
            # supersonic branch: f increasing from 0 at M=1.
            pass
        if supersonic:
            if f(mid) < 0: lo = mid
            else: hi = mid
        else:
            if f(mid) > 0: lo = mid
            else: hi = mid
    return 0.5*(lo+hi)

def read_vtk(path):
    lines = open(path).read().split('\n')
    nx, ny = None, None
    for l in lines:
        if l.startswith('DIMENSIONS'):
            nx, ny, _ = [int(v) for v in l.split()[1:]]
    def field(name):
        i = lines.index(f"SCALARS {name} double 1")
        vals = []; j = i + 2
        while len(vals) < nx*ny and j < len(lines):
            vals += [float(v) for v in lines[j].split()]; j += 1
        return vals[:nx*ny], nx, ny
    return field

def run(tag, vtk):
    field = read_vtk(vtk)
    m, nx, ny = field('mach')
    p, _, _ = field('pressure')
    rho, _, _ = field('density')
    idx = lambda i, j: j*nx + i
    print(f"=== {tag} ===")
    # centerline mach vs 1d area-mach every 10 cells
    print("  x      M_cl    M_1d(sub)  M_1d(sup)")
    x1 = 10.0
    for i in range(0, nx, 16):
        x = (i + 0.5) * x1 / nx
        ar = area(x) / area(5.0)
        msub = mach_from_area_ratio(ar, False)
        msup = mach_from_area_ratio(ar, True)
        print(f"  {x:5.2f}  {m[idx(i,0)]:6.3f}  {msub:8.3f}  {msup:8.3f}")
    # mass flux at inlet/throat/exit on the centerline plane: sum over j of rho*u*A_cell
    dx = x1 / nx
    y1 = None
    # grid y: read from the run dir; the case uses y1=1.204 (uniform). recover from spacing:
    for l in open(vtk).read().split('\n'):
        if l.startswith('SPACING'):
            dx, dy, _ = [float(v) for v in l.split()[1:]]
            y1 = dy * ny
    def mdot(i):
        s = 0.0
        for j in range(ny):
            k = idx(i, j)
            r_in = (j) * dy; r_out = (j+1) * dy
            s += rho[k] * m[k] * math.pi * (r_out**2 - r_in**2)
        return s
    # mach here = |u|/a; u = m*a — a = sqrt(g p/rho). recompute u properly:
    def u_from(k):
        a = math.sqrt(g * max(p[k],1e-12) / max(rho[k],1e-12))
        return m[k] * a
    def mdot2(i):
        s = 0.0
        for j in range(ny):
            k = idx(i, j)
            if rho[k] < 1e-9:
                continue
            r_in = j * dy; r_out = (j+1) * dy
            s += rho[k] * u_from(k) * math.pi * (r_out**2 - r_in**2)
        return s
    mi, mt, me = mdot2(0), mdot2(nx//2), mdot2(nx-1)
    print(f"  mass flux: inlet {mi:.5f}  throat {mt:.5f}  exit {me:.5f}")
    print(f"  spread (max-min)/mean: {(max(mi,mt,me)-min(mi,mt,me))/((mi+mt+me)/3)*100:.2f}%")
    # exit plane state + thrust coefficient (only meaningful for the choked runs)
    pe = p[idx(nx-1, 0)]
    me_cl = m[idx(nx-1, 0)]
    pt = 1.0
    print(f"  exit centerline: M={me_cl:.4f} p/pt={pe/pt:.4f}")

if __name__ == "__main__":
    for tag, vtk in [("p16", "p16.vtk"), ("p89", "p89.vtk"), ("p75", "p75.vtk")]:
        try:
            run(tag, vtk)
        except FileNotFoundError:
            print(f"=== {tag}: {vtk} missing, skipped")

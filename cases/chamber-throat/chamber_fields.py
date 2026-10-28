"""the chamber case's field readers, without the plotting.

the figures and the gate both need the case file, the vtks and the
watcher's history parsed the same way, and the gate runs in the
continuous integration where matplotlib is not installed. so the readers
live here, on the standard library alone, and the figures import them.

readers:
    read_case      the case file's mesh, wall and initial state
    read_vtk       one vtk snapshot into flat field arrays
    columns        the duct's per-column averages
    throat_column  the column holding the wall's narrowest point
    read_history   the watcher's throat history
"""
import ast
import math
import os


def read_case(path):
    """the case's mesh, the wall's polygon and the chamber's state."""
    text = open(path).read()
    sections = {}
    current = None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            current = s[1:-1]
            sections[current] = {}
            continue
        if current is None or "=" not in s:
            continue
        k, v = s.split("=", 1)
        sections[current][k.strip()] = v.split("#")[0].strip()
    mesh = sections["mesh"]
    body = sections["body"]
    ic = sections["initial_condition"]
    return {
        "nx": int(mesh["nx"]),
        "ny": int(mesh["ny"]),
        "x0": float(mesh.get("x0", 0.0)),
        "x1": float(mesh["x1"]),
        "y0": float(mesh.get("y0", 0.0)),
        "y1": float(mesh.get("y1", 0.1)),
        "verts": ast.literal_eval(body["verts"]),
        "p0": float(ic["p0"]),
        "t0": float(ic["t0"]),
    }

def wall_radius(verts, x):
    """the wall's radius at x: the lower of the polygon's crossings."""
    best = float("inf")
    n = len(verts)
    for i in range(n):
        ax, ay = verts[i]
        bx, by = verts[(i + 1) % n]
        if (ax - x) * (bx - x) > 0.0:
            continue
        t = 0.0 if abs(bx - ax) < 1e-12 else (x - ax) / (bx - ax)
        if 0.0 <= t <= 1.0:
            best = min(best, ay + t * (by - ay))
    return best

def signed_distance(verts, px, py):
    """the signed distance to the polygon: negative inside."""
    best = float("inf")
    inside = False
    n = len(verts)
    for i in range(n):
        ax, ay = verts[i]
        bx, by = verts[(i + 1) % n]
        ex, ey = bx - ax, by - ay
        l2 = max(ex * ex + ey * ey, 1e-30)
        t = max(0.0, min(1.0, ((px - ax) * ex + (py - ay) * ey) / l2))
        best = min(best, math.hypot(px - ax - t * ex, py - ay - t * ey))
        if (ay > py) != (by > py) and px < ax + (py - ay) * ex / ey:
            inside = not inside
    return -best if inside else best

def solid_mask(case):
    """the solver's rule: a cell within 0.45 min(dx,dy) of the wall is
    solid, so the body fattens by half a cell."""
    nx, ny = case["nx"], case["ny"]
    dx = (case["x1"] - case["x0"]) / nx
    dy = (case["y1"] - case["y0"]) / ny
    h = min(dx, dy)
    mask = []
    for j in range(ny):
        for i in range(nx):
            x = case["x0"] + (i + 0.5) * dx
            y = case["y0"] + (j + 0.5) * dy
            mask.append(signed_distance(case["verts"], x, y) <= 0.45 * h)
    return mask

def read_vtk(path):
    """the snapshot's fields and its grid."""
    lines = open(path).read().splitlines()
    nx = ny = 0
    origin = (0.0, 0.0)
    spacing = (1.0, 1.0)
    for line in lines:
        parts = line.split()
        if parts[:1] == ["DIMENSIONS"]:
            nx, ny = int(parts[1]), int(parts[2])
        elif parts[:1] == ["ORIGIN"]:
            origin = (float(parts[1]), float(parts[2]))
        elif parts[:1] == ["SPACING"]:
            spacing = (float(parts[1]), float(parts[2]))
    fields = {}
    i = 0
    count = nx * ny
    while i < len(lines):
        if lines[i].startswith("SCALARS"):
            name = lines[i].split()[1]
            values = []
            j = i + 1
            while j < len(lines) and len(values) < count:
                for tok in lines[j].split():
                    try:
                        values.append(float(tok))
                    except ValueError:
                        pass
                j += 1
            if len(values) == count:
                fields[name] = values
            i = j
        else:
            i += 1
    return fields, nx, ny, origin, spacing

def columns(case, fields, nx, ny, origin, spacing):
    """the per-column state: the area-weighted means and the mass flow
    over the fluid rows only."""
    mask = solid_mask(case)
    dx, dy = spacing
    rho = fields["density"]
    pressure = fields["pressure"]
    temperature = fields["temperature"]
    mach = fields["mach"]
    velocity_x = fields.get("velocity_x")
    out = {}
    for i in range(nx):
        area = mdot = psum = tsum = msum = 0.0
        for j in range(ny):
            if mask[j * nx + i]:
                continue
            k = j * nx + i
            rm = origin[1] + j * dy
            rp = rm + dy
            a = math.pi * (rp * rp - rm * rm)
            area += a
            if velocity_x is not None:
                mdot += rho[k] * velocity_x[k] * a
            psum += pressure[k] * a
            tsum += temperature[k] * a
            msum += mach[k] * a
        if area <= 0.0:
            continue
        out[i] = {
            "area": area,
            "m_dot": mdot,
            "p": psum / area,
            "T": tsum / area,
            "mach": msum / area,
        }
    return out

def throat_column(case, cols):
    """the cell holding the wall's narrowest point.

    the grid's throat is a staircase plateau several cells wide, and the
    flow only reaches mach one at the plateau's downstream edge, which
    is where the wall's narrowest vertex sits.
    """
    narrowest = float("inf")
    x_throat = 0.0
    for x, r in case["verts"]:
        if r < narrowest:
            narrowest = r
            x_throat = x
    dx = (case["x1"] - case["x0"]) / case["nx"]
    return min(case["nx"] - 1, int((x_throat - case["x0"]) / dx)), narrowest

def read_history(path):
    """the watcher's history: the numeric columns, in order."""
    rows = []
    with open(path) as f:
        header = f.readline().strip().split(",")
        for line in f:
            if not line.strip():
                continue
            values = line.split(",")
            row = {}
            for key, value in zip(header, values):
                try:
                    row[key] = float(value)
                except ValueError:
                    continue
            rows.append(row)
    return rows

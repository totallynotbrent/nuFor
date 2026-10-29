"""the chamber case's field readers, without the plotting.

the figures and the gate both need the case file, the vtks and the
watcher's history parsed the same way, and the gate runs in the
continuous integration where matplotlib is not installed. so the readers
live here, on the standard library alone, and the figures import them.

readers:
    read_case      the case file's mesh, wall and initial state
    read_vtk       one vtk snapshot into flat field arrays
    columns        the duct's per-column averages
    station_state  a column's mixed mach, flux mach and velocity
    sonic_point    where the marched flow reaches mach one
    range_state    the pooled state over a range of columns
    flux_composition  the gas the mass flux carries
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

_MASK_CACHE = {}


def solid_mask(case):
    """the solver's rule: a cell is solid when its CENTRE lies inside the
    body.

    the solver dropped the older half-cell fattening (body.rs caches the
    mask on `d <= 0.0`), because swallowing the near-wall fluid thickened
    the passage and moved the effective throat and exit area. this reader
    had kept the old rule and was therefore measuring a different
    throat than the solver marched, which is where the stale throat-area
    figure came from.

    the mask is cached on the wall's shape, because the readers ask for
    it once per column and the distance work is the slow part. the rule
    is part of the key, so a future change of rule cannot be served from
    a stale entry.
    """
    key = (
        "centre-inside",
        case["nx"],
        case["ny"],
        case["x0"],
        case["x1"],
        case["y0"],
        case["y1"],
        tuple(tuple(vertex) for vertex in case["verts"]),
    )
    cached = _MASK_CACHE.get(key)
    if cached is not None:
        return cached
    nx, ny = case["nx"], case["ny"]
    dx = (case["x1"] - case["x0"]) / nx
    dy = (case["y1"] - case["y0"]) / ny
    mask = []
    for j in range(ny):
        for i in range(nx):
            x = case["x0"] + (i + 0.5) * dx
            y = case["y0"] + (j + 0.5) * dy
            mask.append(signed_distance(case["verts"], x, y) <= 0.0)
    _MASK_CACHE[key] = mask
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


def station_state(case, fields, nx, ny, origin, spacing, col):
    """the column's station state.

    two mach numbers come out of a column. the mixed mach squares the
    speeds and divides by the reciprocals, and that is the parameter
    the transonic similarity rules use: the area-weighted mean mach is
    no good near the throat, where the wall rows crawl while the core
    runs supersonic. the flux mach is the quasi-one-dimensional one,
    the mass flux's own mach, and that is the one the energy budget
    wants. the station's velocity comes with it.
    """
    mask = solid_mask(case)
    dy = spacing[1]
    rho = fields["density"]
    pressure = fields["pressure"]
    velocity_x = fields.get("velocity_x")
    mach = fields["mach"]
    squared = reciprocal = area = 0.0
    p_sum = t_sum = m_dot = 0.0
    flux = flux_mach_num = flux_vel = 0.0
    for j in range(ny):
        if mask[j * nx + col]:
            continue
        k = j * nx + col
        rm = origin[1] + j * dy
        rp = rm + dy
        a = math.pi * (rp * rp - rm * rm)
        m = max(mach[k], 1.0e-9)
        area += a
        squared += a * m * m
        reciprocal += a / m
        p_sum += a * pressure[k]
        t_sum += a * fields["temperature"][k]
        if velocity_x is not None:
            u = velocity_x[k]
            m_dot += rho[k] * u * a
            weight = rho[k] * u * a
            flux += weight
            flux_mach_num += weight * mach[k]
            flux_vel += weight * u
    if area <= 0.0:
        return None
    return {
        "area": area,
        "p": p_sum / area,
        "T": t_sum / area,
        "m_dot": m_dot,
        "mixed_mach": math.sqrt(squared / reciprocal),
        "flux_mach": (flux_mach_num / flux) if flux > 0.0 else 0.0,
        "velocity": (flux_vel / flux) if flux > 0.0 else 0.0,
    }


def sonic_point(case, fields, nx, ny, origin, spacing):
    """where the marched flow chokes.

    the choke is where the FLUX mach crosses one. the flux mach is the
    mass-weighted mean mach, the same quantity the mass flow balances
    on, and it rises monotonically through the throat, so it gives one
    crossing per run.

    the mixed mach cannot be used for this. the staircase mixes the
    streams inside a cell, so at the throat the mixed mach swings by
    tenths between neighbouring columns (0.53, 0.53, 0.48, 0.32, 0.90,
    0.81, 0.96, 1.09 on 220x80). its first crossing therefore lands on a
    different column on every grid, and the 220x80 record was reading a
    column PAST the throat's single-column plateau.

    the crossing is interpolated between the two straddling columns, and
    it is searched over the whole converging section, because on the
    coarse grids the stream passes mach one slightly AHEAD of the
    geometric minimum: the staircase admits the acceleration early. that
    early crossing is the artefact the throat's geometry has to answer
    for, so it is reported where it happens rather than hidden.
    """
    dx, _ = spacing
    states = {}
    for col in range(nx):
        state = station_state(case, fields, nx, ny, origin, spacing, col)
        if state is not None:
            states[col] = state
    for col in sorted(states):
        if col == 0:
            continue
        left = states[col - 1]
        right = states[col]
        if left["flux_mach"] < 1.0 <= right["flux_mach"]:
            frac = (1.0 - left["flux_mach"]) / max(right["flux_mach"] - left["flux_mach"], 1.0e-12)
            return {
                "col": col,
                "frac": frac,
                "x": case["x0"] + (col - 1 + 0.5 + frac) * dx,
                "p": left["p"] + frac * (right["p"] - left["p"]),
                "T": left["T"] + frac * (right["T"] - left["T"]),
                "m_dot": left["m_dot"] + frac * (right["m_dot"] - left["m_dot"]),
                "area": left["area"] + frac * (right["area"] - left["area"]),
                "basis": "flux_mach",
            }
    return None


def range_state(case, fields, nx, ny, origin, spacing, columns):
    """the pooled state over a range of columns, cell by cell."""
    mask = solid_mask(case)
    dy = spacing[1]
    squared = reciprocal = area = 0.0
    p_sum = t_sum = m_dot = 0.0
    flux = flux_mach_num = flux_vel = 0.0
    for col in columns:
        for j in range(ny):
            if mask[j * nx + col]:
                continue
            k = j * nx + col
            rm = origin[1] + j * dy
            rp = rm + dy
            a = math.pi * (rp * rp - rm * rm)
            m = max(fields["mach"][k], 1.0e-9)
            area += a
            squared += a * m * m
            reciprocal += a / m
            p_sum += a * fields["pressure"][k]
            t_sum += a * fields["temperature"][k]
            if fields.get("velocity_x") is not None:
                u = fields["velocity_x"][k]
                m_dot += fields["density"][k] * u * a
                weight = fields["density"][k] * u * a
                flux += weight
                flux_mach_num += weight * fields["mach"][k]
                flux_vel += weight * u
    if area <= 0.0:
        return None
    return {
        "area": area,
        "p": p_sum / area,
        "T": t_sum / area,
        "m_dot": m_dot,
        "mixed_mach": math.sqrt(squared / reciprocal),
        "flux_mach": (flux_mach_num / flux) if flux > 0.0 else 0.0,
        "velocity": (flux_vel / flux) if flux > 0.0 else 0.0,
    }



def composition(case, fields, nx, ny, origin, spacing, columns, names):
    """the area-weighted mass fractions over a range of columns."""
    mask = solid_mask(case)
    dy = spacing[1]
    sums = {name: 0.0 for name in names}
    area = 0.0
    for col in columns:
        for j in range(ny):
            if mask[j * nx + col]:
                continue
            k = j * nx + col
            rm = origin[1] + j * dy
            rp = rm + dy
            a = math.pi * (rp * rp - rm * rm)
            for name in names:
                key = f"y_{name}"
                if key in fields:
                    sums[name] += a * fields[key][k]
            area += a
    if area <= 0.0:
        return None
    return {name: value / area for name, value in sums.items()}

def flux_composition(case, fields, nx, ny, origin, spacing, columns, names):
    """the mass-flux weighted mass fractions over a range of columns.

    the chamber is a recirculating cavity: the gas on the way to the
    throat is not the gas the chamber holds everywhere. the streamtube
    that feeds the nozzle is the one the mass flux picks out, so the
    composition that decides what the throat chokes at is weighted by
    the mass flux and not by the area.
    """
    mask = solid_mask(case)
    dy = spacing[1]
    sums = {name: 0.0 for name in names}
    flux = 0.0
    for col in columns:
        for j in range(ny):
            if mask[j * nx + col]:
                continue
            k = j * nx + col
            rm = origin[1] + j * dy
            rp = rm + dy
            a = math.pi * (rp * rp - rm * rm)
            weight = fields["density"][k] * fields.get("velocity_x", [0.0] * (nx * ny))[k] * a
            flux += weight
            for name in names:
                key = f"y_{name}"
                if key in fields:
                    sums[name] += weight * fields[key][k]
    if flux <= 0.0:
        return None
    return {name: value / flux for name, value in sums.items()}

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

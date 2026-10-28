"""follow a chamber-throat march and keep its interval snapshots.

the cli overwrites its snapshot file every interval, so a run that is
watched in flight is the only way to keep the settling history. this
copies each new snapshot aside and writes the throat's state beside it.

    python3 watch.py <run-dir> <case.toml> [seconds]

the run dir is where the cli writes: the snapshot file
`<name>.vtk` and the case's own log. the watcher leaves three things
there:

    throat-history.csv   step,t,throat_col,throat_area,throat_mdot,
                         throat_p,throat_T,throat_M,chamber_mdot,
                         chamber_p,chamber_T,chamber_M
    vtk-XXXXXX.vtk       one snapshot per interval
    watch.log            the watcher's own notes

the figures script reads exactly these.
"""
import ast
import csv
import json
import os
import re
import sys
import time

HEADER = [
    "step",
    "t",
    "throat_col",
    "throat_area",
    "throat_mdot",
    "throat_p",
    "throat_T",
    "throat_M",
    "chamber_mdot",
    "chamber_p",
    "chamber_T",
    "chamber_M",
]


def read_case(path):
    """the case's mesh, the wall's polygon and its name."""
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
    return {
        "name": sections["metadata"]["name"].strip('"'),
        "nx": int(mesh["nx"]),
        "ny": int(mesh["ny"]),
        "x0": float(mesh.get("x0", 0.0)),
        "x1": float(mesh["x1"]),
        "y0": float(mesh.get("y0", 0.0)),
        "y1": float(mesh.get("y1", 0.1)),
        "verts": ast.literal_eval(body["verts"]),
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
    """the signed distance to the wall polygon: negative inside."""
    best = float("inf")
    inside = False
    n = len(verts)
    for i in range(n):
        ax, ay = verts[i]
        bx, by = verts[(i + 1) % n]
        ex, ey = bx - ax, by - ay
        l2 = max(ex * ex + ey * ey, 1e-30)
        t = max(0.0, min(1.0, ((px - ax) * ex + (py - ay) * ey) / l2))
        best = min(best, ((px - ax - t * ex) ** 2 + (py - ay - t * ey) ** 2) ** 0.5)
        if (ay > py) != (by > py):
            if px < ax + (py - ay) * ex / ey:
                inside = not inside
    return -best if inside else best


def read_fields(path):
    """the vtk's scalar fields and its grid."""
    lines = open(path).read().splitlines()
    nx = ny = 0
    for line in lines:
        if line.startswith("DIMENSIONS"):
            nx, ny = [int(v) for v in line.split()[1:3]]
    out = {}
    i = 0
    while i < len(lines):
        if lines[i].startswith("SCALARS"):
            name = lines[i].split()[1]
            vals = []
            j = i + 1
            while j < len(lines) and len(vals) < nx * ny:
                for tok in lines[j].split():
                    try:
                        vals.append(float(tok))
                    except ValueError:
                        pass
                j += 1
            if len(vals) == nx * ny:
                out[name] = vals
            i = j
        else:
            i += 1
    return out, nx, ny


def throat_column(case):
    """the cell holding the wall's narrowest point."""
    narrowest = float("inf")
    x_throat = 0.0
    for x, r in case["verts"]:
        if r < narrowest:
            narrowest = r
            x_throat = x
    dx = (case["x1"] - case["x0"]) / case["nx"]
    return min(case["nx"] - 1, int((x_throat - case["x0"]) / dx))


def station(case, col, fields, fluid):
    """a column's area-weighted state and its mass flow."""
    nx, ny = case["nx"], case["ny"]
    dy = (case["y1"] - case["y0"]) / ny
    rho = fields["density"]
    pressure = fields["pressure"]
    temperature = fields["temperature"]
    mach = fields["mach"]
    velocity_x = fields["velocity_x"]
    area = m_dot = p_sum = t_sum = m_sum = 0.0
    for j in range(ny):
        if not fluid[j]:
            continue
        k = j * nx + col
        rm = case["y0"] + j * dy
        ring = 3.141592653589793 * ((rm + dy) ** 2 - rm**2)
        area += ring
        m_dot += rho[k] * velocity_x[k] * ring
        p_sum += pressure[k] * ring
        t_sum += temperature[k] * ring
        m_sum += mach[k] * ring
    return {
        "area": area,
        "m_dot": m_dot,
        "p": p_sum / area,
        "T": t_sum / area,
        "M": m_sum / area,
    }


def march_clock(path, fallback):
    """the march's step and time from the cli's progress file."""
    try:
        with open(path) as fh:
            data = json.load(fh)
        return int(data.get("step", fallback)), float(data.get("t", 0.0))
    except (OSError, ValueError, TypeError):
        return fallback, 0.0


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    run_dir = sys.argv[1]
    case = read_case(sys.argv[2])
    budget = float(sys.argv[3]) if len(sys.argv) > 3 else 7200.0

    nx, ny = case["nx"], case["ny"]
    dx = (case["x1"] - case["x0"]) / nx
    dy = (case["y1"] - case["y0"]) / ny
    erode = 0.45 * min(dx, dy)
    xs = [case["x0"] + (i + 0.5) * dx for i in range(nx)]
    fluid = []
    for j in range(ny):
        y = case["y0"] + (j + 0.5) * dy
        fluid.append([signed_distance(case["verts"], x, y) > erode for x in xs])

    ti = throat_column(case)
    ci = max(1, ti // 2)
    source = os.path.join(run_dir, f"{case['name']}.vtk")
    progress = os.path.join(run_dir, "progress.json")
    history = os.path.join(run_dir, "throat-history.csv")
    log = open(os.path.join(run_dir, "watch.log"), "a")
    seen = set()
    if os.path.exists(history):
        with open(history) as fh:
            seen = {row["vtk"] for row in csv.DictReader(fh) if row.get("vtk")}
    deadline = time.time() + budget
    print(f"watching {source} for the throat at column {ti}", file=log, flush=True)
    while time.time() < deadline:
        if not os.path.exists(source):
            time.sleep(5.0)
            continue
        stamp = os.path.getmtime(source)
        if stamp in seen:
            time.sleep(5.0)
            continue
        try:
            fields, gnx, gny = read_fields(source)
        except (OSError, ValueError):
            time.sleep(5.0)
            continue
        if gnx != nx or gny != ny:
            time.sleep(5.0)
            continue
        seen.add(stamp)
        step, t = march_clock(progress, len(seen))
        throat = station(case, ti, fields, [fluid[j][ti] for j in range(ny)])
        chamber = station(case, ci, fields, [fluid[j][ci] for j in range(ny)])
        new = os.path.join(run_dir, f"vtk-{step:06d}.vtk")
        with open(new, "w") as out:
            out.write(open(source).read())
        write_header = not os.path.exists(history) or os.path.getsize(history) == 0
        with open(history, "a", newline="") as fh:
            writer = csv.writer(fh, lineterminator="\n")
            if write_header:
                writer.writerow(HEADER + ["vtk"])
            writer.writerow(
                [
                    step,
                    f"{t:.3e}",
                    ti,
                    f"{throat['area']:.9e}",
                    f"{throat['m_dot']:.6f}",
                    f"{throat['p']:.6f}",
                    f"{throat['T']:.4f}",
                    f"{throat['M']:.6f}",
                    f"{chamber['m_dot']:.6f}",
                    f"{chamber['p']:.6f}",
                    f"{chamber['T']:.4f}",
                    f"{chamber['M']:.6f}",
                    os.path.basename(new),
                ]
            )
        print(
            f"step {step}: throat M {throat['M']:.4f} p {throat['p']/1e6:.4f} MPa, "
            f"chamber M {chamber['M']:.4f} p {chamber['p']/1e6:.4f} MPa",
            file=log,
            flush=True,
        )
        time.sleep(5.0)
    log.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())

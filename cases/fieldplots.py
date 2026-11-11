"""house-style field figures for every published 2d result raster.

reads each case's own legacy vtk snapshot next to its case.toml and plots
each field on the white-background template shared with the chemistry
figures: filled contours on labelled metre axes with equal aspect, the body
drawn as a white polygon with a dark outline, and a slim vertical colorbar
(short side, pad=0.02 shrink=0.55 aspect=22) carrying units. schlieren is
the viewer's own definition (crates/nufor-cli/src/casework.rs): cellwise
ln(1+|grad rho|), gradient on cell centres, percentiled 2/98 over fluid
cells.

case-specific faithful details:

- the fay case's vtk header spacing is fictional boilerplate for a clustered
  file mesh, so the renderer rebuilds axes from mesh_faces.txt instead.
- the curvilinear o-grid snapshot is STRUCTURED_GRID: cell centroids ride in
  the POINTS block and scalars align with them j*nx+i.
- snapshots that omit a temperature scalar (the o-grid capsule exports
  rho/p/energy/mach only) get T from the checked-in CEA closure inverted on
  the energy axis: bisect log10 e so p(rho, e) matches the exported p, then
  read T off the same table. the nondimensional nozzle lanes use T = p/rho.
- solid cells are masked with the solver's own analytic bodies (sphere /
  sphere-cone / polygon containment), matching crates/nufor-core body.rs.

usage:
    env -u PYTHONPATH /usr/bin/python3 cases/fieldplots.py <out-dir> [--only STEM ...]
"""
import math
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.path import Path as MplPath

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figtemplate import figure, save  # noqa: E402

CASES = os.path.join(HERE, "cases")
DOC_IMG = os.path.join(HERE, "docs", "img")
CEP_TBL = os.path.join(HERE, "crates", "nufor-core", "data", "cea_eqair_table.tsv")

CMAPS = {"rho": "viridis", "p": "plasma", "t": "inferno",
         "mach": "cividis", "schlieren": "gray", "yn": "magma"}
LABELS = {"rho": ("density", "kg/m^3"), "p": ("pressure", "Pa"),
          "t": ("temperature", "K"), "mach": ("Mach number", "-"),
          "yn": ("atomic nitrogen mass fraction", "-"),
          "schlieren": ("numerical schlieren", "ln(1+|grad rho|)")}
PHYS_UNITS = {"rho": ("density", "-"), "p": ("pressure", "-"),
              "t": ("temperature = p/rho", "-"), "mach": ("Mach number", "-"),
              "yn": ("atomic nitrogen mass fraction", "-"),
              "schlieren": LABELS["schlieren"]}
BODY_FACE = "white"
BODY_EDGE = "#2B2B2B"
ORDER = ("mach", "rho", "p", "t", "yn", "schlieren")
MARGIN_X = 0.62
CB_GAP = 0.45
CB_WIDTH = 0.16
TITLE_SPACE = 0.95
LABEL_SPACE = 0.65


def layout(sx, sy):
    """figure size and axes rectangles in inches, honouring equal aspect.

    the data box keeps the domain's own proportions (inches, not axes
    fractions) so the equal-aspect squeeze never shrinks the field, and
    the colorbar is anchored to that box at 55% of its height.
    """
    dw = min(9.0, 6.5 * sx / max(sy, 1e-9))
    dh = dw * sy / max(sx, 1e-9)
    fw = dw + MARGIN_X + CB_GAP + CB_WIDTH + 0.4
    fh = dh + TITLE_SPACE + LABEL_SPACE
    axes = (MARGIN_X / fw, LABEL_SPACE / fh, dw / fw, dh / fh)
    cb_h = 0.55 * dh / fh
    cax = ((MARGIN_X + dw + CB_GAP) / fw,
           (LABEL_SPACE + 0.225 * dh) / fh, CB_WIDTH / fw, cb_h)
    return fw, fh, axes, cax

# ------------------------------------------------------------------ vtk io
def _blocks(text, pos=0):
    out = {}
    for m in re.finditer(r"^SCALARS (\w+) double 1\s*\nLOOKUP_TABLE \w+\s*\n",
                         text, re.M):
        if m.start() < pos:
            continue
        end = text.find("\nSCALARS", m.end())
        stop = text.find("\nCELL_DATA", m.end())
        if end == -1 or (0 < stop < end):
            end = stop if stop != -1 else len(text)
        out[m.group(1)] = np.fromstring(text[m.end():end], sep=" ")
    return out


def read_snapshot(path):
    with open(path, encoding="utf-8", errors="ignore") as fh:
        text = fh.read()
    nx, ny, _ = (int(v) for v in re.search(r"DIMENSIONS (\d+) (\d+) (\d+)",
                                           text).groups())
    ds = re.search(r"DATASET (\w+)", text).group(1)
    if ds == "STRUCTURED_GRID":
        pm = re.search(r"POINTS (\d+) double", text)
        npts = int(pm.group(1))
        start = text.index("\n", pm.end()) + 1
        end = text.find("\nPOINT_DATA", start)
        if end == -1:
            end = text.find("\nSCALARS", start)
        xyz = np.fromstring(text[start:end], sep=" ")[: npts * 3].reshape(npts, 3)
        raw = _blocks(text, pos=start)
        fields = {k: v[:npts].reshape(ny, nx) for k, v in raw.items()}
        return dict(kind="curved", nx=nx, ny=ny,
                    x=xyz[:, 0].reshape(ny, nx), y=xyz[:, 1].reshape(ny, nx)), fields
    ox, oy, _ = (float(v) for v in re.search(r"ORIGIN\s+(\S+)\s+(\S+)\s+(\S+)",
                                             text).groups())
    dx, dy, _ = (float(v) for v in re.search(r"SPACING\s+(\S+)\s+(\S+)\s+(\S+)",
                                             text).groups())
    raw = _blocks(text)
    fields = {k: v[: nx * ny].reshape(ny, nx) for k, v in raw.items()}
    return dict(kind="rect", nx=nx, ny=ny, x0=ox, y0=oy, dx=dx, dy=dy), fields


def fay_faces(case_dir, nx, ny):
    vals = [float(l) for l in open(os.path.join(case_dir, "mesh_faces.txt"))
            if l.strip()]
    return np.array(vals[: nx + 1]), np.array(vals[nx + 1: nx + ny + 2])


# ---------------------------------------------------- eqair closure (rust twin)
NR, NE = 40, 60
LR_G = np.array([-9.231584, -8.911658, -8.591733, -8.271807, -7.951882,
                 -7.631956, -7.312031, -6.992105, -6.672179, -6.352254,
                 -6.032328, -5.712403, -5.392477, -5.072552, -4.752626,
                 -4.432701, -4.112775, -3.792850, -3.472924, -3.152998,
                 -2.833073, -2.513147, -2.193222, -1.873296, -1.553371,
                 -1.233445, -0.913520, -0.593594, -0.273668, 0.046257,
                 0.366183, 0.686108, 1.006034, 1.325959, 1.645885,
                 1.965810, 2.285736, 2.605662, 2.925587, 3.245513])
LE_G = np.array([5.155448, 5.201190, 5.246931, 5.292672, 5.338413, 5.384154,
                 5.429895, 5.475636, 5.521377, 5.567119, 5.612860, 5.658601,
                 5.704342, 5.750083, 5.795824, 5.841565, 5.887306, 5.933048,
                 5.978789, 6.024530, 6.070271, 6.116012, 6.161753, 6.207494,
                 6.253235, 6.298977, 6.344718, 6.390459, 6.436200, 6.481941,
                 6.527682, 6.573423, 6.619164, 6.664906, 6.710647, 6.756388,
                 6.802129, 6.847870, 6.893611, 6.939352, 6.985094, 7.030835,
                 7.076576, 7.122317, 7.168058, 7.213799, 7.259540, 7.305281,
                 7.351023, 7.396764, 7.442505, 7.488246, 7.533987, 7.579728,
                 7.625469, 7.671210, 7.716952, 7.762693, 7.808434, 7.854175])

_ln = open(CEP_TBL).read().splitlines()
TT = np.zeros((NR, NE))
PP = np.zeros((NR, NE))
_i = 3
for _r in range(NR):
    for _c in range(NE):
        _a, _b = _ln[_i].split()[:2]
        TT[_r, _c] = float(_a)
        PP[_r, _c] = float(_b)
        _i += 1


def p_t_at(rho, e):
    r = np.clip(np.log10(np.maximum(rho, 1e-12)), LR_G[0], LR_G[-1])
    ee = np.clip(np.log10(np.maximum(e, 1e-3)), LE_G[0], LE_G[-1])
    i = np.clip(np.searchsorted(LR_G, r, side="right") - 1, 0, NR - 2)
    j = np.clip(np.searchsorted(LE_G, ee, side="right") - 1, 0, NE - 2)
    fi = np.where(LR_G[i + 1] > LR_G[i], (r - LR_G[i]) / (LR_G[i + 1] - LR_G[i]), 0.0)
    fj = np.where(LE_G[j + 1] > LE_G[j], (ee - LE_G[j]) / (LE_G[j + 1] - LE_G[j]), 0.0)
    t = (TT[i, j] * (1 - fi) * (1 - fj) + TT[i + 1, j] * fi * (1 - fj)
         + TT[i, j + 1] * (1 - fi) * fj + TT[i + 1, j + 1] * fi * fj)
    p = (PP[i, j] * (1 - fi) * (1 - fj) + PP[i + 1, j] * fi * (1 - fj)
         + PP[i, j + 1] * (1 - fi) * fj + PP[i + 1, j + 1] * fi * fj)
    return t, p


def eqair_temperature(rho, p):
    """bisect log10 e so closure p(rho, e) matches p; return closure T."""
    r = np.maximum(rho, 1e-12)
    tgt = np.maximum(p, 1e-3)
    lo = np.full_like(r, LE_G[0])
    hi = np.full_like(r, LE_G[-1])
    em = np.power(10.0, 0.5 * (lo + hi))
    for _ in range(48):
        _, pm = p_t_at(r, em)
        lt = pm < tgt
        lo = np.where(lt, np.log10(em), lo)
        hi = np.where(lt, hi, np.log10(em))
        em = np.power(10.0, 0.5 * (lo + hi))
    t, _ = p_t_at(r, em)
    return t


# ------------------------------------------------------------- toml bits
def toml_section(case_dir, name, fname="case.toml"):
    path = os.path.join(case_dir, fname)
    if not os.path.exists(path):
        return ""
    text = open(path, encoding="utf-8").read()
    m = re.search(r"(?ms)^\[" + re.escape(name) + r"\]\s*\n(.*?)(?=^\[|\Z)", text)
    return m.group(1) if m else ""


def kv(section, key, cast=float):
    m = re.search(r"^\s*" + re.escape(key) + r'\s*=\s*"?([^"\n]+)"?',
                  section, re.M)
    return cast(m.group(1).strip()) if m else None


def body_spec(case_dir, fname="case.toml"):
    sec = toml_section(case_dir, "body", fname)
    kind = kv(sec, "type", str) or ""
    if kind == "sphere":
        return dict(kind="sphere", r=kv(sec, "r"),
                    cx=kv(sec, "cx"), cy=kv(sec, "cy"))
    if kind == "sphere_cone":
        return dict(kind="sphere_cone", rn=kv(sec, "rn"),
                    delta=math.radians(kv(sec, "delta_deg")), rb=kv(sec, "rb"),
                    xc=kv(sec, "xc"))
    if kind == "polygon":
        vm = re.search(r"verts\s*=\s*\[(.*?)\]\s*$", sec, re.S | re.M)
        if vm:
            pairs = re.findall(r"\[\s*(-?[\d.eE+-]+)\s*,\s*(-?[\d.eE+-]+)\s*\]",
                               vm.group(1))
            if len(pairs) >= 3:
                return dict(kind="polygon", verts=np.array(
                    [[float(a), float(b)] for a, b in pairs]))
    return dict(kind="none")


# ----------------------------------------------------------- body shapes
def outline_of(spec, n=360):
    k = spec["kind"]
    if k == "sphere":
        th = np.linspace(0, 2 * math.pi, n)
        return np.column_stack([spec["cx"] + spec["r"] * np.cos(th),
                                spec["cy"] + spec["r"] * np.sin(th)])
    if k == "sphere_cone":
        rn, dl, rb, xc = spec["rn"], spec["delta"], spec["rb"], spec["xc"]
        tan = math.tan(dl)
        xt = xc + rn * math.sin(dl)
        yt = rn * math.cos(dl)
        xb = xt + (rb - yt) / tan
        xs = np.linspace(xc - rn, xt, n // 2)
        cap = np.column_stack([xs, np.sqrt(np.maximum(rn**2 - (xs - xc)**2, 0.0))])
        xf = np.linspace(xt, xb, n // 2)
        flank = np.column_stack([xf, yt + (xf - xt) * tan])
        top = np.vstack([cap, flank])
        bottom = top[::-1] * [1.0, -1.0]
        return np.vstack([top, [xb, 0.0], bottom])
    if k == "polygon":
        return spec["verts"]
    return None


def inside_body(spec, X, Y):
    k = spec["kind"]
    if k == "sphere":
        return np.hypot(X - spec["cx"], Y - spec["cy"]) < spec["r"]
    if k == "sphere_cone":
        rn, dl, rb, xc = spec["rn"], spec["delta"], spec["rb"], spec["xc"]
        tan = math.tan(dl)
        sin, cos = math.sin(dl), math.cos(dl)
        xt = xc + rn * sin
        yt = rn * cos
        xb = xt + (rb - yt) / tan
        r_hat = np.hypot(X - xc, Y)
        theta = np.arctan2(np.abs(Y), X - xc)
        cone_half = math.pi / 2 - dl
        nose = (theta > cone_half) & (r_hat < rn)
        flank_zone = (X >= xt) & (X <= xb)
        y_line = yt + (X - xt) * tan
        flank = flank_zone & (np.abs(Y) < y_line) & (np.abs(Y) <= rb)
        cyl = (X > xb) & (np.abs(Y) < rb)
        return nose | flank | cyl
    if k == "polygon":
        pts = np.column_stack([X.ravel(), Y.ravel()])
        return MplPath(spec["verts"]).contains_points(pts).reshape(X.shape)
    return np.zeros(X.shape, bool)


# --------------------------------------------------------------- fields
def grad_mag(field, dx, dy):
    gy, gx = np.gradient(field, dy, dx)
    gx[:, 0] = gx[:, 1] = 0.0
    gy[0, :] = gy[-1, :] = 0.0
    return np.log1p(np.hypot(gx, gy))


def prepare_fields(geom, fields, case_dir, rect_coords):
    phys = toml_section(case_dir, "physics")
    eos = kv(phys, "eos", str) or ""
    out = {}
    alias = {"density": "rho", "pressure": "p", "temperature": "t",
             "y_N": "yn"}
    for name, arr in fields.items():
        out[alias.get(name, name)] = np.array(arr, dtype=float)
    for arr in out.values():
        arr[~np.isfinite(arr)] = np.nan
        arr[np.abs(arr) > 1e13] = np.nan
    synth_t = False
    if "t" not in out and "rho" in out and "p" in out:
        if eos == "eqair":
            out["t"] = eqair_temperature(out["rho"], out["p"])
        else:
            r = kv(phys, "gas_constant")
            r = r if (r is not None and r > 1.0) else 1.0
            den = np.where(np.abs(out["rho"]) > 1e-12, out["rho"], np.nan)
            out["t"] = out["p"] / (den * r)
        synth_t = True
    if "rho" in out and "schlieren" not in out and rect_coords is not None:
        xe, ye = rect_coords
        out["schlieren"] = schlieren_cells(out["rho"], xe, ye)
    return out


def schlieren_cells(rho, xe, ye):
    """central-difference |grad rho| on possibly non-uniform cell edges."""
    gx = np.zeros_like(rho)
    gy = np.zeros_like(rho)
    xcm = 0.5 * (xe[:-1] + xe[1:])
    ycm = 0.5 * (ye[:-1] + ye[1:])
    if rho.shape[1] > 2:
        num = rho[:, 2:] - rho[:, :-2]
        den = (xcm[2:] - xcm[:-2])[None, :]
        gx[:, 1:-1] = num / den
    if rho.shape[0] > 2:
        num = rho[2:, :] - rho[:-2, :]
        den = (ycm[2:] - ycm[:-2])[:, None]
        gy[1:-1, :] = num / den
    return np.log1p(np.hypot(gx, gy))


# --------------------------------------------------------------- panels
def proj_schlieren(XC, YC, RHO, nx=480, ny=320):
    import matplotlib.tri as mtri
    xr = np.asarray(XC).ravel()
    yr = np.asarray(YC).ravel()
    rv = np.asarray(RHO).ravel()
    gi_x = np.linspace(xr.min(), xr.max(), nx)
    gi_y = np.linspace(yr.min(), yr.max(), ny)
    GX, GY = np.meshgrid(gi_x, gi_y)
    tri = mtri.Triangulation(xr, yr)
    interp = mtri.LinearTriInterpolator(tri, rv)
    RG = interp(GX, GY).data
    gy, gx = np.gradient(RG, gi_y, gi_x)
    gx[:, 0] = gx[:, 1] = 0.0
    gy[0, :] = gy[-1, :] = 0.0
    mag = np.log1p(np.hypot(gx, gy))
    back = mtri.LinearTriInterpolator(mtri.Triangulation(GX.ravel(), GY.ravel()),
                                      mag.ravel())
    return back(XC, YC)


def draw_rect(ax, cax, xe, ye, Z, fname, unitless):
    lab, un = PHYS_UNITS[fname] if unitless else LABELS[fname]
    Xe, Ye = np.meshgrid(xe, ye)
    zz = np.ma.masked_invalid(Z)
    fin = zz.compressed()
    lo, hi = ((np.percentile(fin, 2), np.percentile(fin, 98))
              if fin.size else (0.0, 1.0))
    if hi - lo < 1e-12:
        hi = lo + 1.0
    im = ax.pcolormesh(Xe, Ye, zz, cmap=CMAPS[fname], vmin=lo, vmax=hi,
                       shading="flat", rasterized=True)
    cb = plt.colorbar(im, cax=cax)
    cb.set_label(f"{lab} ({un})", fontsize=10)
    cb.ax.tick_params(labelsize=9)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal", adjustable="box")


def draw_curved(ax, cax, XC, YC, Z, fname, unitless):
    import matplotlib.tri as mtri
    lab, un = PHYS_UNITS[fname] if unitless else LABELS[fname]
    tri = mtri.Triangulation(np.asarray(XC).ravel(), np.asarray(YC).ravel())
    zv = np.asarray(Z).ravel()
    fin = zv[np.isfinite(zv)]
    lo, hi = ((np.percentile(fin, 2), np.percentile(fin, 98))
              if fin.size else (0.0, 1.0))
    if hi - lo < 1e-12:
        hi = lo + 1.0
    cmap = plt.get_cmap(CMAPS[fname]).copy()
    cmap.set_bad(alpha=0.0)
    mp = np.ma.masked_invalid(zv)
    tc = ax.tripcolor(tri, mp, cmap=cmap, vmin=lo, vmax=hi, shading="flat",
                      rasterized=True)
    cb = plt.colorbar(tc, cax=cax)
    cb.set_label(f"{lab} ({un})", fontsize=10)
    cb.ax.tick_params(labelsize=9)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_aspect("equal", adjustable="box")


BASE_NAME = {"rho": "Density", "p": "Pressure", "t": "Temperature",
             "mach": "Mach number", "yn": "Atomic nitrogen mass fraction",
             "schlieren": "Numerical schlieren"}
METAS = {
    "aeroshell-h30": ("aeroshell capsule, trim 22.5 deg", "h = 30.0 km, v = 4274 m/s, equilibrium air"),
    "aeroshell-h34.6": ("aeroshell capsule, trim 22.5 deg", "h = 34.6 km, v = 6691 m/s, equilibrium air"),
    "aeroshell-h40": ("aeroshell capsule, trim 22.5 deg", "h = 40.0 km, v = 8706 m/s, equilibrium air"),
    "aeroshell-h42.5": ("aeroshell capsule, trim 22.5 deg", "h = 42.5 km, v = 9327 m/s, equilibrium air"),
    "capsule-eqair": ("entry capsule Rn = 1.9 m", "v = 6960 m/s, equilibrium air, cartesian 320x160"),
    "capsule-eqair-fine": ("entry capsule Rn = 1.9 m", "v = 6960 m/s, equilibrium air, cartesian 640x320"),
    "capsule-ogrid": ("entry capsule Rn = 1.9 m", "v = 6960 m/s, equilibrium air, body-fitted O-grid"),
    "cyl-fay": ("cylinder R = 1 m", "M = 20 N2 freestream v = 6047 m/s, equilibrium air, viscous"),
    "cyl-park5": ("Mach 20 N2 cylinder R = 1 m", "reacting, park 5-species kinetics"),
    "cyl-park5visc-full": ("Mach 20 N2 cylinder R = 1 m", "park 5-species, fully catalytic 1000 K wall"),
    "cyl-park5visc-none": ("Mach 20 N2 cylinder R = 1 m", "park 5-species, non-catalytic adiabatic wall"),
    "nozzle-p16": ("CDV nozzle cold flow", "reservoir ratio p_exit/p0 = 0.16, perfect gas"),
    "nozzle-p75": ("CDV nozzle cold flow", "reservoir ratio p_exit/p0 = 0.75, perfect gas"),
    "nozzle-p89": ("CDV nozzle cold flow", "reservoir ratio p_exit/p0 = 0.89, perfect gas"),
}

SPECS = [
    ("aeroshell-h30.0-a22.5", "aeroshell-h30.0-a22.5.vtk", "aeroshell-h30", "case.toml"),
    ("aeroshell-h34.6-a22.5", "aeroshell-h34.6-a22.5.vtk", "aeroshell-h34.6", "case.toml"),
    ("aeroshell-h40.0-a22.5", "aeroshell-h40.0-a22.5.vtk", "aeroshell-h40", "case.toml"),
    ("aeroshell-h42.5-a22.5", "aeroshell-h42.5-a22.5.vtk", "aeroshell-h42.5", "case.toml"),
    ("capsule-reentry-eqair", "capsule-reentry-eqair.vtk", "capsule-eqair", "case.toml"),
    ("capsule-reentry-eqair-fine", "capsule-reentry-eqair-fine.vtk", "capsule-eqair-fine", "case.toml"),
    ("capsule-reentry-ogrid", "capsule-reentry-ogrid.vtk", "capsule-ogrid", "case.toml"),
    ("cylinder-stagnation-fay", "cylinder-stagnation-fay.vtk", "cyl-fay", "case.toml"),
    ("cylinder-m20-park5", "cylinder-m20-park5.vtk", "cyl-park5", "case.toml"),
    ("cylinder-m20-park5-visc", "cylinder-m20-park5-visc-full.vtk", "cyl-park5visc-full", "case-full.toml"),
    ("cylinder-m20-park5-visc", "cylinder-m20-park5-visc-none.vtk", "cyl-park5visc-none", "case-none.toml"),
    ("nozzle-cdv", "nozzle-cdv-p16.vtk", "nozzle-p16", "case-p16.toml"),
    ("nozzle-cdv", "nozzle-cdv-p75.vtk", "nozzle-p75", "case-p75.toml"),
    ("nozzle-cdv", "nozzle-cdv-p89.vtk", "nozzle-p89", "case-p89.toml"),
]


def render_one(out_dir, case_rel, vtk_name, stem, cfg):
    case_dir = os.path.join(CASES, case_rel)
    geom, fields = read_snapshot(os.path.join(case_dir, vtk_name))
    spec = body_spec(case_dir, cfg)
    outline = outline_of(spec)
    unitless = "nozzle" in stem

    if geom["kind"] == "rect":
        nx, ny = geom["nx"], geom["ny"]
        if case_rel == "cylinder-stagnation-fay":
            xe, ye = fay_faces(case_dir, nx, ny)
        else:
            xe = geom["x0"] + np.arange(nx + 1) * geom["dx"]
            ye = geom["y0"] + np.arange(ny + 1) * geom["dy"]
        XC, YC = np.meshgrid(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]))
        rect_coords = (xe, ye)
    else:
        XC, YC = geom["x"], geom["y"]
        xe = np.linspace(XC.min(), XC.max(), geom["nx"] + 1)
        ye = np.linspace(YC.min(), YC.max(), geom["ny"] + 1)
        rect_coords = None

    fields = prepare_fields(geom, fields, case_dir, rect_coords)
    if rect_coords is None and "schlieren" not in fields and "rho" in fields:
        fields["schlieren"] = proj_schlieren(XC, YC, fields["rho"])
    solid = inside_body(spec, XC, YC)
    if "mach" in fields:
        solid |= ~np.isfinite(fields["mach"]) | (fields["mach"] == 0.0)

    sx = float(XC.max() - XC.min())
    sy = float(YC.max() - YC.min())
    fw, fh, axes_box, cb_box = layout(sx, sy)
    subj, cond = METAS[stem]

    for fname in ORDER:
        if fname not in fields:
            continue
        Z = np.array(fields[fname], dtype=float)
        Z[solid] = np.nan
        fig, ax = figure(figsize=(fw, fh))
        ax.set_position(axes_box)
        cax = fig.add_axes(cb_box)
        if rect_coords is not None:
            draw_rect(ax, cax, xe, ye, Z, fname, unitless)
        else:
            draw_curved(ax, cax, XC, YC, Z, fname, unitless)
        if outline is not None:
            ax.add_patch(MplPolygon(outline, closed=True, facecolor=BODY_FACE,
                                    edgecolor=BODY_EDGE, lw=1.3, zorder=6))
        ax.set_title(f"{BASE_NAME[fname]}, {subj}, {cond}", fontsize=11)
        ax.autoscale_view()
        save(fig, f"{stem}-{fname}", out_dir=out_dir)


def main(out_dir, only=None):
    os.makedirs(out_dir, exist_ok=True)
    for case_rel, vtk_name, stem, cfg in SPECS:
        if only and not any(stem.startswith(o) for o in only):
            continue
        print("family", stem, flush=True)
        render_one(out_dir, case_rel, vtk_name, stem, cfg)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else DOC_IMG
    only = None
    if "--only" in sys.argv:
        only = sys.argv[sys.argv.index("--only") + 1:]
    main(out, only)

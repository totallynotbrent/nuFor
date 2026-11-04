"""reports the throat's choked state and its geometry.

the sonar point is the choke: the place where the flux mach crosses one
after the last minimum-area column. the mixed mach oscillates around the
throat because the staircase mixes the streams across a cell, so it picks
different columns on different grids; the flux mach rises monotonically
through the throat and is the same quantity the mass flow balances on, so
it is the crossing the choke contract should use.

the geometry comes from the solver's own rule, the same way the marcher
classifies its cells, so the reported throat is the marched throat.

usage: throat_state.py <run-dir> <case.toml>
"""
import glob
import os
import sys

sys.path.insert(0, "/home/ebarleb/Bots/nuFor/cases/chamber-throat")
import chamber_fields as cf  # noqa: E402


def build_mask(case):
    """the solver's mask: a cell is solid when its centre is inside."""
    nx, ny = case["nx"], case["ny"]
    dx = (case["x1"] - case["x0"]) / nx
    dy = (case["y1"] - case["y0"]) / ny
    verts = [tuple(vertex) for vertex in case["verts"]]
    return [
        cf.signed_distance(verts, (i + 0.5) * dx, (j + 0.5) * dy) <= 0.0
        for j in range(ny)
        for i in range(nx)
    ]


def column_states(case, fields, nx, ny, origin, spacing):
    """the per-column state of every column with fluid."""
    states = {}
    for col in range(nx):
        state = cf.station_state(case, fields, nx, ny, origin, spacing, col)
        if state is not None:
            states[col] = state
    return states


def main():
    run_dir = sys.argv[1]
    case_path = sys.argv[2]
    case = cf.read_case(case_path)
    nx = case["nx"]
    dx = (case["x1"] - case["x0"]) / nx

    vtks = sorted(glob.glob(os.path.join(run_dir, "vtk-*.vtk")))
    fields, fnx, fny, origin, spacing = cf.read_vtk(vtks[-1])

    # the tail of the run, averaged the way the distiller averages it
    tails = []
    for path in vtks[-8:]:
        tf, tfx, tfy, torigin, tspacing = cf.read_vtk(path)
        state = cf.sonic_point(case, tf, tfx, tfy, torigin, tspacing)
        if state is not None:
            tails.append(state)
    mixed_mean = sum(item["p"] for item in tails) / len(tails) / 1e6 if tails else float("nan")

    states = column_states(case, fields, fnx, fny, origin, spacing)
    finite = {col: state["area"] for col, state in states.items() if state["area"] > 0.0}
    minimum = min(finite.values())
    minima = sorted(col for col, area in finite.items() if abs(area - minimum) < 1e-15)
    last = minima[-1]

    # the choke: the flux mach crossing one. on the coarse grids the
    # crossing lands INSIDE the minimum-area plateau, because the
    # staircase lets the stream accelerate past mach one before the
    # geometric minimum. so the search runs over the whole converging
    # section and takes the FIRST crossing, which is the choke.
    crossing = None
    for col in sorted(states):
        if col == 0:
            continue
        left = states[col - 1]["flux_mach"]
        right = states[col]["flux_mach"]
        if left < 1.0 <= right:
            frac = (1.0 - left) / max(right - left, 1.0e-12)
            crossing = (col - 1, col, frac)
            break
    if crossing is None:
        # already supersonic at the first column: report the first column
        crossing = (0, 1, 0.0)

    print("case %s grid %dx%d" % (os.path.basename(run_dir), nx, case["ny"]))
    print(
        "  throat_area_m2 %.6e  min_columns %s  last_min %d  throat_x_m %.6f"
        % (minimum, minima, last, (last + 0.5) * dx)
    )
    if crossing is None:
        print("  the flow never reaches flux mach one after the throat")
        return
    left_col, right_col, frac = crossing
    left = states[left_col]
    right = states[right_col]
    pressure = left["p"] + frac * (right["p"] - left["p"])
    print(
        "  choke between cols %d and %d (frac %.4f) x_m %.6f  (geometric throat at col %d)"
        % (left_col, right_col, frac, (left_col + 0.5 + frac) * dx, last)
    )
    print(
        "  pressure_at_choke_mpa %.4f   mixed_mach_detector_mpa %.4f"
        % (pressure / 1e6, mixed_mean)
    )
    print("  col   area        p_MPa     Mmix     Mflux")
    for col in range(max(0, last - 6), min(nx, last + 8)):
        state = states.get(col)
        if state is None:
            continue
        print(
            "  %3d %.4e %9.4f %7.4f %8.5f"
            % (col, state["area"], state["p"] / 1e6, state["mixed_mach"], state["flux_mach"])
        )


if __name__ == "__main__":
    main()

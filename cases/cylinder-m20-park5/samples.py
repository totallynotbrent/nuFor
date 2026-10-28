# distil the marched field of the dissociation cylinder into the samples
# the gate reads, so the published numbers can be audited from a clean
# checkout without marching the case again.
#
# run from the case directory:
#   nufor run case.toml         # writes cylinder-m20-park5.vtk
#   python3 samples.py          # writes samples.csv
#
# the vtk itself is a run product and stays out of the repository. this
# script writes the slices the gate measures: the freestream slab, the
# stagnation line, and the atomic nitrogen map of the whole field.
import csv
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "case.toml")
VTK = os.path.join(HERE, "cylinder-m20-park5.vtk")
SAMPLES = os.path.join(HERE, "samples.csv")

# the slice widths and the section names the gate expects.
FREESTREAM_WIDTH_DIVISOR = 30
SECTION_INLET = "inlet"
SECTION_STAGNATION = "stagnation"
SECTION_FIELD = "field_n"

# the scalar fields the samples carry, in the order they are written.
SLICE_FIELDS = (
    "density",
    "pressure",
    "temperature",
    "mach",
    "y_N",
    "y_N2",
    "y_O2",
    "y_O",
    "y_NO",
)


def read_case(path):
    """the case's own mesh extent."""
    text = open(path).read()

    def grab(key):
        m = re.search(r"(?m)^\s*" + key + r"\s*=\s*([-\d.eE+]+)", text)
        return float(m.group(1)) if m else None

    return {"nx": grab("nx"), "ny": grab("ny"), "x0": grab("x0"),
            "x1": grab("x1"), "r": grab("r"), "cx": grab("cx")}


def read_fields(path):
    """every scalar field in the vtk, row major with x fastest."""
    text = open(path).read()
    dim = re.search(r"DIMENSIONS (\d+) (\d+)", text)
    if not dim:
        raise ValueError(f"{path}: no DIMENSIONS block")
    nx, ny = int(dim.group(1)), int(dim.group(2))
    n = nx * ny
    fields = {}
    for name, body in re.findall(
        r"SCALARS (\S+) double 1\nLOOKUP_TABLE default\n((?:[-\d.eE+ ]+\n)+)", text
    ):
        vals = [float(tok) for tok in body.split()]
        if len(vals) >= n:
            fields[name] = vals[:n]
    if not fields:
        raise ValueError(f"{path}: no scalar fields")
    return fields, nx, ny


def fmt(value):
    """six significant digits, enough to audit and small enough to ship."""
    return f"{value:.6g}"


def main():
    if not os.path.exists(VTK):
        print(f"MISS: no {os.path.basename(VTK)}, run the case first")
        return 1
    case = read_case(CASE)
    fields, nx, ny = read_fields(VTK)
    missing = [f for f in SLICE_FIELDS if f not in fields]
    if missing:
        print(f"FAIL: the field carries no {missing}")
        return 1
    width = max(2, nx // FREESTREAM_WIDTH_DIVISOR)

    with open(SAMPLES, "w", newline="") as fh:
        out = csv.writer(fh)
        out.writerow(["# park5 samples"])
        out.writerow(["# section", SECTION_INLET])
        out.writerow(["i"] + list(SLICE_FIELDS))
        for i in range(width):
            out.writerow([i] + [fmt(fields[f][i]) for f in SLICE_FIELDS])

        out.writerow(["# section", SECTION_STAGNATION])
        out.writerow(["i"] + list(SLICE_FIELDS))
        for i in range(nx):
            out.writerow([i] + [fmt(fields[f][i]) for f in SLICE_FIELDS])

        out.writerow(["# section", SECTION_FIELD])
        out.writerow(["j", "y_N"])
        for j in range(ny):
            row = ",".join(fmt(fields["y_N"][j * nx + i]) for i in range(nx))
            out.writerow([j, row])

    print(
        f"samples: {width} freestream cells, {nx} stagnation cells, "
        f"{nx}x{ny} atomic nitrogen map"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

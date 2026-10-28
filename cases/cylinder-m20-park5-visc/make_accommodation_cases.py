# the catalytic-wall accommodation sweep: one case per accommodation
# coefficient, generated from the shipped non-catalytic case so the
# sweep changes nothing but the wall condition.
#
#   python3 make_accommodation_cases.py [out_dir]
#
# each generated case carries `catalytic_wall = "partial"` and a
# `catalytic_efficiency` in [0, 1]; the ends of the ladder (0 and 1)
# reproduce the shipped none and full cases.
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "case-none.toml")

ACCOMMODATION = (0.0, 0.25, 0.5, 0.75, 1.0)


def case_text(eta):
    label = f"{eta:g}"
    if eta == 0.0:
        condition = "none"
    elif eta == 1.0:
        condition = "full"
    else:
        condition = "partial"
    text = open(BASE).read()
    text = text.replace(
        'name = "cylinder-m20-park5-visc-none"',
        f'name = "cylinder-m20-park5-visc-eta-{label}"',
    ).replace(
        'description = "Mach 20 N2 cylinder, park5 viscous lane, '
        'non-catalytic isothermal 1000 K wall: the q_w bracket lower bound"',
        f'description = "Mach 20 N2 cylinder, park5 viscous lane, '
        f'isothermal 1000 K wall, accommodation {label}"',
    ).replace('catalytic_wall = "none"', f'catalytic_wall = "{condition}"')
    if condition == "partial":
        text = text.replace(
            f'catalytic_wall = "{condition}"',
            f'catalytic_wall = "{condition}"\ncatalytic_efficiency = {label}',
        )
    return text


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "accommodation")
    os.makedirs(out_dir, exist_ok=True)
    for eta in ACCOMMODATION:
        label = f"{eta:g}"
        path = os.path.join(out_dir, f"case-eta-{label}.toml")
        with open(path, "w") as fh:
            fh.write(case_text(eta))
        print(path)


if __name__ == "__main__":
    main()

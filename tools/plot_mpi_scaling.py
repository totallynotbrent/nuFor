#!/usr/bin/env python3
"""draw the mpi scaling figures in docs/numerics/parallel from a bench json.

    tools/plot_mpi_scaling.py SCALING_JSON --out docs/numerics/parallel

produces mpi-scaling.png (speedup against ideal) and, when the json rows
come from the viscous lane, mpi-viscous-timing.png (bars of march wall
time). the figures are derived artifacts: the page's table cites the
same json, so text, table, and picture cannot drift apart.

run it with an interpreter that has matplotlib: on fedora that is
/usr/bin/python3, not the default python3 on some setups.

    /usr/bin/python3 tools/plot_mpi_scaling.py SCALING_JSON \
        --out docs/numerics/parallel
"""

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("json")
    parser.add_argument("--out", default="docs/numerics/parallel")
    args = parser.parse_args()

    with open(args.json) as handle:
        payload = json.load(handle)
    rows = payload["rows"]
    ranks = [r["ranks"] for r in rows]
    best = [r["best_s"] for r in rows]
    speedup = [r["speedup"] for r in rows]

    os.makedirs(args.out, exist_ok=True)
    if payload.get("lane") == "viscous":
        fig, ax = plt.subplots(figsize=(5.2, 3.4), dpi=150)
        ax.bar([str(r) for r in ranks], best, color="#4878cf", width=0.6)
        for idx, (tick, val) in enumerate(zip(ax.patches, best)):
            rect = ax.patches[idx]
            ax.text(
                rect.get_x() + rect.get_width() / 2,
                val,
                f"{val:.1f}s",
                ha="center",
                va="bottom",
                fontsize=8,
            )
        ax.set_xlabel("ranks (np)")
        ax.set_ylabel("march wall time (s)")
        ax.set_title("laminar viscous lane, 320x176 eqair march")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "mpi-viscous-timing.png"))
        plt.close(fig)
        print("wrote mpi-viscous-timing.png")
        return

    fig, ax = plt.subplots(figsize=(5.6, 3.6), dpi=150)
    ax.plot(ranks, ranks, "--", color="#999999", label="ideal linear")
    ax.plot(ranks, speedup, "o-", color="#4878cf", label="measured speedup")
    for x, y in zip(ranks, speedup):
        ax.annotate(
            f"{y:.2f}x",
            (x, y),
            textcoords="offset points",
            xytext=(0, 9),
            ha="center",
            fontsize=8,
        )
    ax.set_xticks(ranks)
    ax.set_ylim(0, max(ranks) + 0.6)
    ax.set_xlabel("ranks (np)")
    ax.set_ylabel("speedup vs np=1")
    ax.set_title("strong scaling, 320x176 eqair euler march")
    ax.legend(fontsize=8, loc="upper left", framealpha=0.95)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "mpi-scaling.png"))
    plt.close(fig)
    print("wrote mpi-scaling.png")


if __name__ == "__main__":
    main()

"""shared figure template for all nuFor docs/research figures.

usage: from a case script, `from figtemplate import figure; fig, ax =
figure(nrows=1)` or `fig, axes = figure(nrows=3)`. style the axes,
save with figtemplate.save(fig, name) which writes
../../docs/research/<name>.png (doc-relative) at 150 dpi.
"""
import os
import inspect
from cycler import cycler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# the palette every published figure shares
ORANGE = "#D35400"
PURPLE = "#7D74B5"
GREEN = "#16A085"
RED = "#C0392B"
GRAY = "#A9A9A9"

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.grid": False,
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "legend.facecolor": "white",
    "axes.prop_cycle": cycler(color=[ORANGE, PURPLE, GREEN, RED, GRAY]),
})


def figure(nrows=1, ncols=1, **kwargs):
    """new figure on the house template. size defaults to a 2:1 canvas."""
    if "figsize" not in kwargs:
        kwargs["figsize"] = (7 * ncols, 3.6 * nrows)
    return plt.subplots(nrows, ncols, **kwargs)


def save(fig, name, out_dir=None):
    """write docs/research/<name>.png relative to the CALLER's case dir.

    interactive callers (stdin, -c strings) have no real file, so they
    pass out_dir explicitly.
    """
    if out_dir is None:
        caller = os.path.dirname(os.path.abspath(inspect.stack()[1].filename))
        out_dir = os.path.normpath(os.path.join(caller, "..", "..", "docs", "research"))
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, name + ".png")
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", path)
    return path

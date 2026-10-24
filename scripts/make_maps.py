"""Map figure for the article: observed 1985 and 1999, and the 1999 maps simulated
by lulcc and by the replication (CLUE-S, deterministic level).

Usage: python scripts/make_maps.py GOLDENS_DIR [--out results/fig_maps.png]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from re_lulcc import clues, data  # noqa: E402

# Colours of the original article's Fig. 3 (palegreen, midnightblue, indianred1).
COLOURS = ["#98fb98", "#191970", "#ff6a6a"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("goldens", type=Path)
    ap.add_argument("--out", type=Path, default=Path("results/fig_maps.png"))
    args = ap.parse_args()
    D = args.goldens / "deterministic"

    pie = data.load_pie(args.goldens / "data")
    gp = pd.read_csv(D / "probabilities.csv.gz")
    gd = pd.read_csv(D / "demand.csv")[data.LABELS].to_numpy()
    prob = {}
    for lab in data.LABELS:
        a = np.zeros(pie.shape)
        a.ravel()[gp["cell"].to_numpy() - 1] = gp[lab].to_numpy()
        prob[lab] = a
    b = data.make_backend(pie, pie.lu[0])
    clues.run_clues(b, prob, gd)
    ours = data.to_codes(b, pie.mask)
    with rasterio.open(D / "clues_maps.tif") as r:
        theirs = r.read(15).astype(np.int16)
    theirs[~pie.mask] = 0
    n_diff = int(((ours != theirs) & pie.mask).sum())

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    cmap = ListedColormap(COLOURS)
    panels = [
        (pie.lu[0], "Observed 1985"),
        (pie.lu[14], "Observed 1999"),
        (theirs, "CLUE-S 1999, lulcc"),
        (ours, "CLUE-S 1999, replication"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(10, 2.9))
    for ax, (m, title) in zip(axes, panels):
        shown = np.ma.masked_where(m == 0, m)
        ax.imshow(shown, cmap=cmap, vmin=1, vmax=3, interpolation="nearest")
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.legend(handles=[Patch(color=c, label=l) for c, l in zip(COLOURS, data.LABELS)],
               loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=250)
    print(f"cells differing between the two simulated maps: {n_diff}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

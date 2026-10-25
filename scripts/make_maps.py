# SPDX-License-Identifier: GPL-3.0-or-later
"""Map figure for the article: observed 1985 and 1999, and the 1999 maps simulated
by lulcc and by the replication with CLUE-S (deterministic level: no perturbation).

The replication's map is produced by its whole chain: its own logistic regression
(statsmodels), its own demand and DisSLUCC's published CLUE-S allocation. From the
reference only the training cells are taken, because R's random sample cannot be
drawn again in Python. lulcc's map is read from the reference (band 15 = 1999).

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

from re_lulcc import clues, data, demand, suitability  # noqa: E402

# Colours of the original article's Fig. 3 (palegreen, midnightblue, indianred1).
COLOURS = ["#98fb98", "#191970", "#ff6a6a"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("goldens", type=Path)
    ap.add_argument("--out", type=Path, default=Path("results/fig_maps.png"))
    args = ap.parse_args()
    D = args.goldens / "deterministic"

    pie = data.load_pie(args.goldens / "data")
    cells, mask = pie.cells, pie.mask
    train = pd.read_csv(D / "train_cells.csv")["cell"].to_numpy() - 1   # the only input from lulcc
    coefs = suitability.fit(pie, train)
    p = suitability.probabilities(pie, coefs, cells)
    prob = {}
    for j, lab in enumerate(data.LABELS):
        a = np.zeros(pie.shape)
        a.ravel()[cells] = p[:, j]
        prob[lab] = a
    totals = np.array([[(pie.lu[t][mask] == c).sum() for c in data.CODES] for t in (0, 6, 14)])
    dmd = demand.interpolated_demand([0, 6, 14], totals, np.arange(15), int(mask.sum()))
    b = data.make_backend(pie, pie.lu[0])
    clues.run_clues(b, prob, dmd)
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

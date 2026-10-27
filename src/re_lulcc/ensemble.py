# SPDX-License-Identifier: GPL-3.0-or-later
"""One run of the paper's experiment with our own randomness (one seed).

Mirrors demo/gmd-paper.R: 10% stratified training partition of the 1985
map, glm per class, AUCs (Figs. 5-6), CLUE-S with jitter 0.0002 and the
stochastic ordered allocation over t = 1..14, and the figure of merit of
forest -> built at t14 at eight resolutions (Fig. 8), with its agreement
budget (Fig. 7).

The partition follows caret::createDataPartition as lulcc calls it: the
class codes are numeric, so caret cuts them at their quantiles
(probs 0, .25, .5, .75, 1; unique breaks; lowest included) and samples
ceiling(10%) of each group. For the Plum Island map this gives two groups,
{forest, built} and {other}.
"""
from __future__ import annotations

import time

import numpy as np
from scipy.stats import rankdata

from . import clues, data, demand, ordered, suitability, threemap

FACTORS = [2 ** k for k in range(1, 9)]
CLUES = dict(elasticity=0.2, scale=1e-6, max_iteration=1000, max_difference=50, jitter=0.0002)
ORDER = [2, 1, 3]  # built, forest, other


def partition(y: np.ndarray, p: float, rng: np.random.Generator) -> np.ndarray:
    """Indices (into y) of the training set, as caret::createDataPartition for numeric y."""
    breaks = np.unique(np.quantile(y, [0, 0.25, 0.5, 0.75, 1.0]))
    group = np.digitize(y, breaks[1:-1], right=True)  # (b0, b1], first bin closed below
    train = []
    for g in np.unique(group):
        idx = np.flatnonzero(group == g)
        n = 1 if len(idx) == 1 else int(np.ceil(len(idx) * p))
        train.append(rng.choice(idx, size=n, replace=False))
    return np.sort(np.concatenate(train))


def auc(score: np.ndarray, label: np.ndarray) -> float:
    """Area under the ROC curve (Mann-Whitney, ties counted as 1/2), as ROCR."""
    label = label.astype(bool)
    n1, n0 = label.sum(), (~label).sum()
    r = rankdata(score)
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def run_seed(pie: data.Pie, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    cells = pie.cells
    lu0 = pie.lu[0].ravel()[cells]

    tr = partition(lu0.astype(float), 0.1, rng)
    test = np.setdiff1d(np.arange(len(cells)), tr)
    coefs = suitability.fit(pie, cells[tr])
    prob = suitability.probabilities(pie, coefs, cells)

    out: dict = {"seed": seed}
    for j, (lab, code) in enumerate(zip(data.LABELS, data.CODES)):
        out[f"auc_glm_{lab}"] = auc(prob[test, j], lu0[test] == code)
    gain = lu0 != 2
    lu6 = pie.lu[6].ravel()[cells]
    out["auc_builtgain"] = auc(prob[gain, 1], lu6[gain] == 2)

    totals = np.array([[(pie.lu[t][pie.mask] == c).sum() for c in data.CODES] for t in (0, 6, 14)])
    dmd = demand.interpolated_demand([0, 6, 14], totals, np.arange(15), int(pie.mask.sum()))

    full = {}
    for j, lab in enumerate(data.LABELS):
        a = np.zeros(pie.shape)
        a.ravel()[cells] = prob[:, j]
        full[lab] = a

    t = time.perf_counter()
    b = data.make_backend(pie, pie.lu[0])
    iters = clues.run_clues(b, full, dmd, rng=rng, **CLUES)
    clues_final = data.to_codes(b, pie.mask)
    out["seconds_clues"] = time.perf_counter() - t
    out["clues_iterations_max"] = int(max(iters))

    t = time.perf_counter()
    om = ordered.run_ordered(prob, lu0.astype(np.int16), dmd, data.CODES, ORDER, stochastic=True, rng=rng)
    ordered_final = np.zeros(pie.shape, dtype=np.int16)
    ordered_final.ravel()[cells] = om[-1]
    out["seconds_ordered"] = time.perf_counter() - t
    out["ordered_unclassified_final"] = int((om[-1] == -1).sum())

    x, x1 = pie.lu[0].astype(np.int16), pie.lu[14].astype(np.int16)
    for name, y1 in (("clues", clues_final), ("ordered", ordered_final)):
        tabs = threemap.tables(x, x1, y1, FACTORS, data.CODES)
        for f, tab in zip(FACTORS, tabs):
            out[f"fom_forest_built_{name}_{f}"] = float(threemap.figure_of_merit(tab, 3)["transition"][0, 1])
        # Agreement budget of forest -> built (Fig. 7 of the paper), as fractions of the
        # study area: a misses, b hits, c wrong hits, d false alarms. Descriptive only:
        # not among the metrics of acceptance.toml.
        for f, tab in zip(FACTORS, tabs):
            ab = threemap.agreement_budget(tab, 3, [0], [1], transition=True)
            for k in "abcd":
                out[f"ab_forest_built_{name}_{f}_{k}"] = float(ab[k])
    return out

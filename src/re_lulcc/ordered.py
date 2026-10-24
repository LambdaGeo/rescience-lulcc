"""Ordered allocation (Sect. 3.5.3 of the paper; Fuchs et al., 2013).

Classes are visited in a fixed order (here built, forest, other). For a class
whose demand rises, the n most suitable cells not yet in it, and not already
settled by an earlier class, are converted. For a class whose demand falls,
its n least suitable cells are released to a temporary class (-1), which
later classes can take. Cells settled by a class are removed from the
competition. In the stochastic version, a cell is a candidate only if its
suitability is >= (rising) or < (falling) a U(0, 1) draw; draws are repeated
until there are enough candidates or 1000 attempts.

Provisional: written from the paper's description; the points the paper
leaves open (which cells count as candidates, the temporary class, the
1000-attempt cap, the order of ties) were resolved against lulcc's
allocate.R and are listed in the README.
"""
from __future__ import annotations

import numpy as np


def ordered_step(prob: np.ndarray, lu0: np.ndarray, demand: np.ndarray, codes: list[int],
                 order: list[int], stochastic: bool, rng: np.random.Generator | None = None) -> np.ndarray:
    """prob: (n_cells, n_classes); lu0: (n_cells,) codes; demand: (n_classes,). Returns new codes."""
    prob = prob.copy()
    area0 = np.array([(lu0 == c).sum() for c in codes])
    if np.abs(demand - area0).sum() == 0:
        return lu0.copy()
    lu1 = lu0.copy()
    for code in order:
        ix = codes.index(code)
        n = int(demand[ix] - (lu1 == code).sum())
        if n == 0:
            prob[lu0 == code, :] = np.nan
        elif n > 0:
            cand = np.flatnonzero(lu1 != code)
            p = prob[cand, ix]
            keep = ~np.isnan(p)
            cand, p = cand[keep], p[keep]
            srt = np.argsort(-p, kind="stable")          # high to low, ties in cell order
            cand, p = cand[srt], p[srt]
            sel = _candidates(p, n, stochastic, rng, rising=True)
            lu1[cand[sel[:n]]] = code
            prob[lu1 == code, :] = np.nan
        else:
            cand = np.flatnonzero(lu0 == code)
            p = prob[cand, ix]
            keep = ~np.isnan(p)
            cand, p = cand[keep], p[keep]
            srt = np.argsort(p, kind="stable")           # low to high
            cand, p = cand[srt], p[srt]
            sel = _candidates(p, -n, stochastic, rng, rising=False)
            lu1[cand[sel[:-n]]] = -1
            prob[lu1 == code, :] = np.nan
    return lu1


def _candidates(p: np.ndarray, n: int, stochastic: bool, rng, rising: bool) -> np.ndarray:
    if not stochastic:
        return np.arange(len(p))
    sel = np.array([], dtype=int)
    for _ in range(1001):
        u = rng.uniform(size=len(p))
        sel = np.flatnonzero(p >= u if rising else p < u)
        if len(sel) >= n:
            break
    return sel


def run_ordered(prob: np.ndarray, lu0: np.ndarray, demand: np.ndarray, codes: list[int],
                order: list[int], stochastic: bool = False,
                rng: np.random.Generator | None = None) -> list[np.ndarray]:
    """All steps; demand rows indexed by t (row 0 = initial). Returns maps for t = 0..T."""
    maps = [lu0.copy()]
    cur = lu0.copy()
    for t in range(1, demand.shape[0]):
        cur = ordered_step(prob, cur, demand[t], codes, order, stochastic, rng)
        maps.append(cur.copy())
    return maps

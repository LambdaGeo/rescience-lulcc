# SPDX-License-Identifier: GPL-3.0-or-later
"""Three-map comparison at multiple resolutions (Pontius et al., 2011).

Compares the initial reference map (x), the final reference map (x1) and the
simulated final map (y1). At each aggregation factor, every coarse cell
holds the fraction of its valid fine cells in each class; the three maps are
then compared fuzzily and summarised in a three-dimensional contingency
table, from which the agreement budget and the figure of merit follow.

Provisional: equations from Pontius et al. (2011). The layout of the table
(a ((n+1)^2, n+1) matrix), the aggregation (blocks from the top-left corner,
partial blocks at the edges kept, weighted by their number of valid cells)
and the treatment of empty denominators follow lulcc's
ThreeMapComparison.R, so that the tables can be compared one to one.

Table indexing (0-based here, 1-based in lulcc): row i + j*(n+1), column k,
for i = class in x, j = class in x1, k = class in y1; index n is "total".
"""
from __future__ import annotations

import numpy as np

# Empty denominators (0/0) are expected where a class is absent; lulcc lets them
# become NaN and drops them in the sums, and so does this module.
_quiet = np.errstate(invalid="ignore", divide="ignore")


def _aggregate(a: np.ndarray, f: int) -> np.ndarray:
    """Block sums ignoring NaN, like raster::aggregate(fun=sum, na.rm=TRUE, expand=TRUE).
    Blocks with no valid cell are NaN."""
    if f == 1:
        return a.copy()
    nr, nc = a.shape
    pr, pc = (-nr) % f, (-nc) % f
    b = np.pad(a, ((0, pr), (0, pc)), constant_values=np.nan)
    b = b.reshape(b.shape[0] // f, f, b.shape[1] // f, f)
    s = np.nansum(b, axis=(1, 3))
    s[np.all(np.isnan(b), axis=(1, 3))] = np.nan
    return s


def tables(x: np.ndarray, x1: np.ndarray, y1: np.ndarray, factors: list[int],
           codes: list[int]) -> list[np.ndarray]:
    """x, x1, y1: integer maps; cells outside the study area are <= 0 or NaN in x."""
    n = len(codes)
    valid = np.isfinite(x.astype(float)) & (x > 0)
    out = []
    for f in factors:
        ones = np.where(valid, 1.0, np.nan)
        w = _aggregate(ones, f).ravel()
        keep = np.isfinite(w)
        w = w[keep]
        Q, R, S = [], [], []
        for c in codes:
            for m, lst in ((x, Q), (x1, R), (y1, S)):
                ind = np.where(valid, (m == c).astype(float), np.nan)
                lst.append(_aggregate(ind, f).ravel()[keep] / w)
        out.append(_table(Q, R, S, w, n))
    return out


def _fmin(*arrs):
    r = arrs[0]
    for a in arrs[1:]:
        r = np.fmin(r, a)
    return r


def _wmean(v: np.ndarray, w: np.ndarray) -> float:
    ok = ~np.isnan(v)
    return float(np.sum(v[ok] * w[ok]) / np.sum(w))


@_quiet
def _table(Q, R, S, w, n) -> np.ndarray:
    eq1 = [_fmin(R[j], S[j]) for j in range(n)]
    eq2 = [_fmin(Q[j], R[j], S[j]) for j in range(n)]

    eq3 = [[None] * n for _ in range(n)]  # eq3[j][i], i != j
    for j in range(n):
        bsum = sum(Q[i] - eq2[i] for i in range(n) if i != j)
        for i in range(n):
            if i != j:
                eq3[j][i] = (eq1[j] - eq2[j]) * (Q[i] - eq2[i]) / bsum

    Qq = []
    for i in range(n):
        parts = [eq3[j][i] if j != i else eq2[i] for j in range(n)]
        T = np.nansum(np.column_stack(parts), axis=1)
        Qq.append(Q[i] - T)
    Rr = [R[j] - eq1[j] for j in range(n)]
    Ss = [S[k] - eq1[k] for k in range(n)]

    a = sum(eq1)
    bsum = sum(Qq[i] * Rr[j] * Ss[k] for j in range(n) for k in range(n) for i in range(n) if j != k)

    tab = np.full(((n + 1) * (n + 1), n + 1), np.nan)
    row = lambda i, j: i + j * (n + 1)  # noqa: E731

    for j in range(n):                                   # rule 1
        tab[row(n, j), j] = _wmean(eq1[j], w)
    for j in range(n):                                   # rule 2
        tab[row(j, j), j] = _wmean(eq2[j], w)
    for j in range(n):                                   # rule 3
        for i in range(n):
            if i != j:
                tab[row(i, j), j] = _wmean(eq3[j][i], w)
    for j in range(n):                                   # rule 4
        for k in range(n):
            if j != k:
                for i in range(n):
                    tab[row(i, j), k] = _wmean((1 - a) * Qq[i] * Rr[j] * Ss[k] / bsum, w)

    for j in range(n):                                   # totals over x1
        for k in range(n):
            tab[j + (n + 1) * n, k] = np.nansum(tab[[row(j, jj) for jj in range(n)], k])
    for j in range(n + 1):                               # totals over x
        for k in range(n):
            tab[row(n, j), k] = np.nansum(tab[[row(i, j) for i in range(n)], k])
    for j in range(n + 1):                               # totals over y1
        for i in range(n + 1):
            tab[row(i, j), n] = np.nansum(tab[row(i, j), :n])
    return tab


@_quiet
def figure_of_merit(tab: np.ndarray, n: int) -> dict:
    """Overall, per class and per transition (Pontius et al., 2011, Eqs. 9-11).
    transition[i, j] is the transition from class i to class j."""
    row = lambda i, j: i + j * (n + 1)  # noqa: E731
    a = np.nansum([tab[row(n, j), j] for j in range(n)])
    b = np.nansum([tab[row(j, j), j] for j in range(n)])
    overall = (a - b) / (1 - b)
    category = np.zeros(n)
    for i in range(n):
        a = np.nansum([tab[row(i, j), j] for j in range(n)])
        bb = np.nansum([tab[row(i, j), n] for j in range(n)])
        c = tab[row(i, i), i]
        category[i] = (a - c) / (bb - c)
    trans = np.zeros((n, n))
    for j in range(n):
        for i in range(n):
            a = tab[row(i, j), j]
            b = tab[row(i, j), n]
            c = np.nansum([tab[row(i, k), j] for k in range(n)])
            trans[i, j] = a / (b + c - a)
    return {"overall": overall, "category": category, "transition": trans}


def agreement_budget(tab: np.ndarray, n: int, from_ix: list[int], to_ix: list[int],
                     transition: bool) -> dict:
    """Sources of agreement and disagreement, as fractions of the study area:
    a change simulated as persistence (misses), b change simulated correctly
    (hits), c change to the wrong class (wrong hits), d persistence simulated
    as change (false alarms), e persistence simulated correctly."""
    row = lambda i, j: i + j * (n + 1)  # noqa: E731
    a = sum(np.nansum([tab[row(i, j), i] for i in from_ix if i != j]) for j in to_ix)
    b = sum(sum(tab[row(i, j), j] for i in from_ix if i != j) for j in to_ix)
    c = 0.0
    for j in to_ix:
        for i in from_ix:
            c += np.nansum([tab[row(i, j), k] for k in range(n) if i != j and i != k and j != k])
    if len(to_ix) == 1:
        (fi,), (ti,) = from_ix, to_ix
        c += sum(tab[row(fi, j), ti] for j in range(n) if j != fi and j != ti)
    d = sum(sum(tab[row(j, j), k] for k in to_ix if k != j) for j in from_ix)
    out = {"a": a, "b": b, "c": c, "d": d}
    if not transition:
        out["e"] = sum(tab[row(i, i), i] for i in from_ix)
    return out

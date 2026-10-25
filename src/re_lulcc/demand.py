# SPDX-License-Identifier: GPL-3.0-or-later
"""Demand by linear interpolation between observed maps (Sect. 3.4 of the paper).

lulcc's approxExtrapDemand() interpolates (or extrapolates) the number of
cells of each class linearly in time, sets negative values to zero, and then
rounds each row to integers that sum to the number of valid cells
(roundSum): floor everything, then give one more cell to the classes with
the largest remainders.
"""
from __future__ import annotations

import numpy as np


def approx_extrap(t_obs: np.ndarray, y_obs: np.ndarray, t_out: np.ndarray) -> np.ndarray:
    """Piecewise-linear interpolation, linear extrapolation from the end segments."""
    t_obs = np.asarray(t_obs, dtype=float)
    y_obs = np.asarray(y_obs, dtype=float)
    t_out = np.asarray(t_out, dtype=float)
    y = np.interp(t_out, t_obs, y_obs)
    lo, hi = t_out < t_obs[0], t_out > t_obs[-1]
    if lo.any():
        slope = (y_obs[1] - y_obs[0]) / (t_obs[1] - t_obs[0])
        y[lo] = y_obs[0] + slope * (t_out[lo] - t_obs[0])
    if hi.any():
        slope = (y_obs[-1] - y_obs[-2]) / (t_obs[-1] - t_obs[-2])
        y[hi] = y_obs[-1] + slope * (t_out[hi] - t_obs[-1])
    return y


def round_sum(x: np.ndarray, ncell: int) -> np.ndarray:
    out = np.empty_like(x, dtype=np.int64)
    for r in range(x.shape[0]):
        y = x[r] / x[r].sum() * ncell
        xint = np.floor(y)
        rem = y - xint
        order = np.argsort(rem, kind="stable")  # ascending, ties keep column order
        tot = int(round(ncell - xint.sum()))
        if tot > 0:
            xint[order[len(y) - tot:]] += 1
        out[r] = xint.astype(np.int64)
    return out


def interpolated_demand(t_obs, totals: np.ndarray, t_out, ncell: int) -> np.ndarray:
    """totals: (n_obs, n_classes) cell counts at t_obs. Returns (len(t_out), n_classes)."""
    cols = [approx_extrap(t_obs, totals[:, j], t_out) for j in range(totals.shape[1])]
    d = np.column_stack(cols)
    d[d < 0] = 0
    return round_sum(d, ncell)

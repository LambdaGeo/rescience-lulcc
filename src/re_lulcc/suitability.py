"""Binary logistic regression per class (Sect. 3.3 of the paper).

One model per class, fitted on the training cells of the 1985 map:
    Built  ~ ef_001 + ef_002 + ef_003
    Forest ~ ef_001 + ef_002
    Other  ~ ef_001 + ef_002
"""
from __future__ import annotations

import numpy as np
import statsmodels.api as sm
from scipy.special import expit

from .data import CODES, LABELS, Pie

FORMULAS: dict[str, list[str]] = {
    "Forest": ["ef_001", "ef_002"],
    "Built": ["ef_001", "ef_002", "ef_003"],
    "Other": ["ef_001", "ef_002"],
}


def fit(pie: Pie, train_cells: np.ndarray, t: int = 0) -> dict[str, dict[str, float]]:
    """train_cells: flat indices (0-based). Returns {class: {"(Intercept)": b0, driver: b}}."""
    lu = pie.lu[t].ravel()[train_cells]
    coefs = {}
    for lab, code in zip(LABELS, CODES):
        cols = FORMULAS[lab]
        X = np.column_stack([pie.drivers[c].ravel()[train_cells] for c in cols])
        X = sm.add_constant(X, has_constant="add")
        y = (lu == code).astype(float)
        res = sm.GLM(y, X, family=sm.families.Binomial()).fit(tol=1e-12, maxiter=100)
        coefs[lab] = {"(Intercept)": float(res.params[0])} | {
            c: float(b) for c, b in zip(cols, res.params[1:])
        }
    return coefs


def probabilities(pie: Pie, coefs: dict[str, dict[str, float]], cells: np.ndarray) -> np.ndarray:
    """(len(cells), 3) suitability, columns in LABELS order, float64."""
    out = np.empty((len(cells), len(LABELS)))
    for j, lab in enumerate(LABELS):
        c = coefs[lab]
        z = np.full(len(cells), c["(Intercept)"])
        for name, beta in c.items():
            if name != "(Intercept)":
                z = z + beta * pie.drivers[name].ravel()[cells]
        out[:, j] = expit(z)
    return out

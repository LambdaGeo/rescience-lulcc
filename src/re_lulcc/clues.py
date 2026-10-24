"""CLUE-S on DisSModel/DisSLUCC (Sect. 3.5.2 of the paper).

The allocation is DisSLUCC's AllocationDClueSLike, unchanged: each cell goes
to the class with the highest suitability plus a per-class iteration term;
the term moves by (demand - allocated) * factor until the largest
difference is within max_difference.

Two pieces are provisional here and are meant to move into DisSLUCC:

* ProbabilityPotential: suitability from a given probability matrix, in
  float64, plus the elasticity of the cell's current class. DisSLUCC's
  PotentialDLogisticRegression computes the same quantity from the
  regression coefficients, but in float32.
* JitteredDClueS: the random perturbation of the paper (`jitter.f`), drawn
  again at every iteration from U(-jitter, jitter) for every cell and
  class, and lulcc's behaviour at the iteration limit (warn and return the
  last allocation instead of raising).
"""
from __future__ import annotations

import warnings

import numpy as np
from disslucc.components.allocation import AllocationDClueSLike
from disslucc.components.demand import DemandPreComputedValues
from dissmodel.core import Environment
from dissmodel.geo import RasterBackend, SyncRasterModel

from .data import LABELS


class ProbabilityPotential(SyncRasterModel):
    """Writes <class>_pot = probability + elasticity * (cell is in that class)."""

    def setup(self, backend, probability: dict[str, np.ndarray], elasticity: dict[str, float]) -> None:
        super().setup(backend)
        self.probability = probability
        self.elasticity = elasticity
        for lab in LABELS:
            self.backend.set(lab + "_pot", np.zeros(self.shape, dtype=np.float64))

    def execute(self) -> None:
        for lab in LABELS:
            current = self.backend.get(lab) == 1
            self.backend.arrays[lab + "_pot"] = (
                self.probability[lab] + np.where(current, self.elasticity[lab], 0.0)
            )


class JitteredDClueS(AllocationDClueSLike):
    """AllocationDClueSLike plus the jitter of lulcc's CLUE-S."""

    def setup(self, backend, *, jitter: float = 0.0, rng: np.random.Generator | None = None,
              on_max_iteration: str = "warn", **kwargs) -> None:
        super().setup(backend, **kwargs)
        self.jitter = float(jitter)
        self.rng = rng if rng is not None else np.random.default_rng()
        self.on_max_iteration = on_max_iteration

    def execute(self) -> None:
        lu_types = self.land_use_types
        n_lu = len(lu_types)
        mask = self._mask()
        regions = np.clip(self.backend.get(self.region_attr).ravel().astype(int) - 1, 0, self._tm.shape[0] - 1)
        start = np.stack([self.backend.get(lu).ravel() for lu in lu_types], axis=1)
        allowed = self._tm[regions, np.argmax(start, axis=1), :]
        pot = np.stack([self.backend.get(lu + "_pot").ravel() for lu in lu_types], axis=1)
        valid = mask.ravel()
        iter_vec = np.zeros(n_lu, dtype=np.float64)

        for n_iter in range(self.max_iteration + 1):
            p = pot
            if self.jitter > 0:
                noise = np.zeros_like(pot)
                noise[valid] = self.rng.uniform(-self.jitter, self.jitter, size=(int(valid.sum()), n_lu))
                p = pot + noise
            scores = np.where(allowed, p + iter_vec[np.newaxis, :], -np.inf)
            best = np.argmax(scores, axis=1)
            for j, lu in enumerate(lu_types):
                new = (best == j).astype(np.float32).reshape(self.shape)
                self.backend.arrays[lu] = np.where(mask, new, self.backend.get(lu))
            diff = self._calc_diff(mask)
            max_diff = max(abs(v) for v in diff.values())
            if max_diff <= self.max_difference:
                self.iterations_per_step.append(n_iter)
                return
            if n_iter >= self.max_iteration:
                msg = (f"step {int(self.env.now())}: no convergence after {self.max_iteration} "
                       f"iterations (max difference {max_diff:.0f} cells)")
                if self.on_max_iteration == "raise":
                    raise RuntimeError(msg)
                warnings.warn(msg, stacklevel=2)
                self.iterations_per_step.append(n_iter)
                return
            for j, lu in enumerate(lu_types):
                iter_vec[j] += diff[lu] * self.factor_iteration


def run_clues(
    backend: RasterBackend,
    probability: dict[str, np.ndarray],
    demand: np.ndarray,
    *,
    elasticity: float = 0.2,
    scale: float = 1e-6,
    max_iteration: int = 1000,
    max_difference: float = 50,
    jitter: float = 0.0,
    rng: np.random.Generator | None = None,
    record=None,
) -> list[int]:
    """Runs t = 1..len(demand)-1 on `backend` in place.

    demand[t] is the number of cells of each class at t (row 0 = initial map);
    lulcc allocates from t = 1, so the environment starts at 1. `record`, if
    given, is called with (t, backend) after every step. Returns the number of
    iterations of each step.
    """
    n = demand.shape[0] - 1
    env = Environment(start_time=1, end_time=n)
    dem = DemandPreComputedValues(
        annual_demand=[list(map(float, row)) for row in demand], land_use_types=LABELS,
    )
    ProbabilityPotential(backend=backend, probability=probability,
                         elasticity={lab: elasticity for lab in LABELS})
    common = dict(
        backend=backend, demand=dem, land_use_types=LABELS,
        transition_matrix=[[[1, 1, 1]] * 3], cell_area=1.0,
        max_difference=max_difference, max_iteration=max_iteration, factor_iteration=scale,
    )
    if jitter == 0 and rng is None:
        alloc = AllocationDClueSLike(**common)   # DisSLUCC 0.5.0 as published
    else:
        alloc = JitteredDClueS(**common, jitter=jitter, rng=rng)
    if record is not None:
        _Recorder(backend=backend, callback=record)
    env.run()
    return alloc.iterations_per_step


class _Recorder(SyncRasterModel):
    def setup(self, backend, callback) -> None:
        super().setup(backend)
        self.callback = callback

    def execute(self) -> None:
        self.callback(int(self.env.now()), self.backend)

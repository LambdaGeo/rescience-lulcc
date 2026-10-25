# SPDX-License-Identifier: GPL-3.0-or-later
"""Plum Island Ecosystems data, read from the GeoTIFFs exported by lulcc-goldens.

Cell order is row-major (C order), the same order as raster's cell numbers,
so a flat index here is lulcc's cell number minus one.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from dissmodel.geo import RasterBackend

LABELS = ["Forest", "Built", "Other"]
CODES = [1, 2, 3]
YEARS = {0: 1985, 6: 1991, 14: 1999}
DRIVERS = ["ef_001", "ef_002", "ef_003"]  # elevation, slope, distance to built (1985)


@dataclass
class Pie:
    lu: dict[int, np.ndarray]        # t -> int8 map, 0 outside the study area
    drivers: dict[str, np.ndarray]   # name -> float64, NaN outside
    mask: np.ndarray                 # bool, True inside the study area
    transform: object
    crs: object

    @property
    def shape(self) -> tuple[int, int]:
        return self.mask.shape

    @property
    def cells(self) -> np.ndarray:
        """Flat indices of valid cells, in lulcc's order."""
        return np.flatnonzero(self.mask.ravel())


def load_pie(data_dir: str | Path) -> Pie:
    data_dir = Path(data_dir)
    lu: dict[int, np.ndarray] = {}
    transform = crs = None
    for t, year in YEARS.items():
        with rasterio.open(data_dir / f"lu_pie_{year}.tif") as r:
            a = r.read(1, masked=True)
            lu[t] = np.where(a.mask, 0, a.filled(0)).astype(np.int8)
            transform, crs = r.transform, r.crs
    drivers = {}
    for name in DRIVERS:
        with rasterio.open(data_dir / f"{name}.tif") as r:
            drivers[name] = r.read(1, masked=True).astype(np.float64).filled(np.nan)
    mask = lu[0] > 0
    return Pie(lu=lu, drivers=drivers, mask=mask, transform=transform, crs=crs)


def one_hot(lu: np.ndarray) -> dict[str, np.ndarray]:
    """Land use as DisSLUCC expects it: one 0/1 array per class."""
    return {lab: (lu == code).astype(np.float32) for lab, code in zip(LABELS, CODES)}


def to_codes(backend: RasterBackend, mask: np.ndarray) -> np.ndarray:
    """One-hot arrays in the backend back to codes 1..3 (0 outside the mask)."""
    out = np.zeros(backend.shape, dtype=np.int8)
    for lab, code in zip(LABELS, CODES):
        out[backend.get(lab) == 1] = code
    out[~mask] = 0
    return out


def make_backend(pie: Pie, lu: np.ndarray) -> RasterBackend:
    b = RasterBackend(shape=pie.shape, transform=pie.transform, crs=pie.crs)
    for name, arr in one_hot(lu).items():
        b.set(name, arr)
    for name, arr in pie.drivers.items():
        b.set(name, np.nan_to_num(arr, nan=0.0))
    b.set("mask", pie.mask.copy())
    return b


def crosstab(a: np.ndarray, b: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """3x3 table, rows = classes in `a`, columns = classes in `b`."""
    out = np.zeros((3, 3), dtype=np.int64)
    for i, ci in enumerate(CODES):
        for j, cj in enumerate(CODES):
            out[i, j] = int(((a == ci) & (b == cj) & mask).sum())
    return out

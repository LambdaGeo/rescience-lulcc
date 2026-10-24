"""Run the replication over seeds and write results/ensemble_python.csv.

Usage: python scripts/run_ensemble.py DATA_DIR [--seeds 30] [--first 1] [--out results/ensemble_python.csv]
DATA_DIR is lulcc-goldens' goldens/data (the pie GeoTIFFs).
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from re_lulcc import data, ensemble  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("--seeds", type=int, default=30)
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--out", type=Path, default=Path("results/ensemble_python.csv"))
    args = ap.parse_args()

    pie = data.load_pie(args.data_dir)
    rows = []
    for s in range(args.first, args.first + args.seeds):
        t = time.perf_counter()
        rows.append(ensemble.run_seed(pie, s))
        print(f"seed {s}: {time.perf_counter() - t:.0f} s", flush=True)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(args.out, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

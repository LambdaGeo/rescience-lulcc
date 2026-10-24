"""Compare the replication with the deterministic golden of lulcc-goldens, stage by stage.

Usage: python scripts/check_deterministic.py GOLDENS_DIR [--out results/deterministic.json]

Stages (each one fed with the golden's own inputs, so a difference points at
that stage only):
  1 data          observed transition tables
  2 glm           coefficients fitted on the golden's training cells
  3 probability   suitability from our coefficients vs the golden's matrix
  4 demand        interpolated, rounded demand
  5 clues         CLUE-S (DisSLUCC, jitter 0) from the golden's probabilities and demand
  6 ordered       ordered allocation (deterministic) from the same inputs
  7 validation    three-map tables and figure of merit from the golden's maps
  8 chain         everything from our own glm onwards
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from re_lulcc import clues, data, demand, ordered, suitability, threemap  # noqa: E402

FACTORS = [2 ** k for k in range(1, 9)]


def read_maps(path: Path, mask: np.ndarray) -> np.ndarray:
    with rasterio.open(path) as r:
        a = r.read().astype(np.int16)
    a[:, ~mask] = 0
    a[a == 255] = 0
    return a


def golden_tables(path: Path) -> dict[int, np.ndarray]:
    df = pd.read_csv(path)
    out = {}
    for f, g in df.groupby("factor"):
        tab = np.full((g["row"].max(), g["col"].max()), np.nan)
        tab[g["row"] - 1, g["col"] - 1] = g["value"]
        out[int(f)] = tab
    return out


def fom_from_golden(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def maps_diff(ours: list[np.ndarray], gold: np.ndarray, mask: np.ndarray) -> list[int]:
    return [int(((o != g) & mask).sum()) for o, g in zip(ours, gold)]


def run_clues_maps(pie, prob_full, dmd):
    b = data.make_backend(pie, pie.lu[0])
    maps = [pie.lu[0].copy()]
    rec = lambda t, be: maps.append(data.to_codes(be, pie.mask))  # noqa: E731
    iters = clues.run_clues(b, prob_full, dmd, record=rec)
    return maps, iters


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("goldens", type=Path)
    ap.add_argument("--out", type=Path, default=Path("results/deterministic.json"))
    args = ap.parse_args()
    G = args.goldens
    D = G / "deterministic"
    report: dict = {}

    pie = data.load_pie(G / "data")
    mask, cells = pie.mask, pie.cells
    n_classes = len(data.LABELS)

    # 1 data
    ct = pd.read_csv(G / "data" / "crosstab.csv")
    ok = True
    for (t0, t1), g in ct.groupby(["t_from", "t_to"]):
        ours = data.crosstab(pie.lu[int(t0)], pie.lu[int(t1)], mask)
        gold = g.pivot(index="from", columns="to", values="cells").loc[data.LABELS, data.LABELS].to_numpy()
        ok &= bool((ours == gold).all())
    report["1_data_crosstab_identical"] = ok

    # 2 glm
    train = pd.read_csv(D / "train_cells.csv")["cell"].to_numpy() - 1
    coefs = suitability.fit(pie, train)
    gc = pd.read_csv(D / "glm_coefficients.csv")
    rel = []
    for _, r in gc.iterrows():
        ours = coefs[r["label"]][r["term"]]
        rel.append(abs(ours - r["estimate"]) / max(abs(r["estimate"]), 1e-300))
    report["2_glm_max_rel_diff"] = float(max(rel))

    # 3 probability
    gp = pd.read_csv(D / "probabilities.csv.gz")
    gcells = gp["cell"].to_numpy() - 1
    assert np.array_equal(gcells, cells), "cell order differs from the golden"
    gprob = gp[data.LABELS].to_numpy()
    rcoefs = {lab: {} for lab in data.LABELS}
    for _, r in gc.iterrows():
        rcoefs[r["label"]][r["term"]] = r["estimate"]
    report["3_prob_from_R_coefs_max_abs_diff"] = float(np.abs(suitability.probabilities(pie, rcoefs, cells) - gprob).max())
    our_prob = suitability.probabilities(pie, coefs, cells)
    report["3_prob_from_our_coefs_max_abs_diff"] = float(np.abs(our_prob - gprob).max())

    # 4 demand
    totals = np.array([[(pie.lu[t][mask] == c).sum() for c in data.CODES] for t in (0, 6, 14)])
    dmd = demand.interpolated_demand([0, 6, 14], totals, np.arange(15), int(mask.sum()))
    gd = pd.read_csv(D / "demand.csv")[data.LABELS].to_numpy()
    report["4_demand_identical"] = bool((dmd == gd).all())

    def full(prob_cells: np.ndarray) -> dict[str, np.ndarray]:
        out = {}
        for j, lab in enumerate(data.LABELS):
            a = np.zeros(pie.shape)
            a.ravel()[cells] = prob_cells[:, j]
            out[lab] = a
        return out

    # 5 CLUE-S
    gclues = read_maps(D / "clues_maps.tif", mask)
    t = time.perf_counter()
    cmaps, iters = run_clues_maps(pie, full(gprob), gd)
    report["5_clues_seconds"] = round(time.perf_counter() - t, 1)
    report["5_clues_cells_differing_per_step"] = maps_diff(cmaps, gclues, mask)
    report["5_clues_iterations_per_step"] = iters

    # 6 ordered
    gord = read_maps(D / "ordered_maps.tif", mask)
    lu0 = pie.lu[0].ravel()[cells].astype(np.int16)
    om = ordered.run_ordered(gprob, lu0, gd, data.CODES, [2, 1, 3], stochastic=False)
    omaps = []
    for m in om:
        a = np.zeros(pie.shape, dtype=np.int16)
        a.ravel()[cells] = m
        omaps.append(a)
    report["6_ordered_cells_differing_per_step"] = maps_diff(omaps, gord, mask)
    report["6_ordered_unclassified_cells_final"] = int((om[-1] == -1).sum())

    # 7 validation, from the golden's own maps
    x, x1 = pie.lu[0].astype(np.int16), pie.lu[14].astype(np.int16)
    for model, gmaps in (("clues", gclues), ("ordered", gord)):
        tabs = threemap.tables(x, x1, gmaps[-1], FACTORS, data.CODES)
        gt = golden_tables(D / f"{model}_tables.csv")
        report[f"7_{model}_tables_max_abs_diff"] = float(max(
            np.nanmax(np.abs(tab - gt[f])) for tab, f in zip(tabs, FACTORS)))
        gf = fom_from_golden(D / f"{model}_fom.csv")
        gfb = gf[(gf.level == "transition") & (gf["from"] == "Forest") & (gf["to"] == "Built")].fom.to_numpy()
        ofb = np.array([threemap.figure_of_merit(tab, n_classes)["transition"][0, 1] for tab in tabs])
        report[f"7_{model}_fom_forest_built_max_abs_diff"] = float(np.abs(ofb - gfb).max())
        report[f"7_{model}_fom_forest_built"] = [round(float(v), 4) for v in ofb]

    # 8 whole chain from our glm
    cm8, _ = run_clues_maps(pie, full(our_prob), dmd)
    report["8_chain_clues_cells_differing_final"] = int(((cm8[-1] != gclues[-1]) & mask).sum())

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2))
    for k, v in report.items():
        print(f"{k:45s} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Apply acceptance.toml to the two ensembles and draw the comparison figure.

Usage: python scripts/compare_ensemble.py R_ENSEMBLE_DIR PY_CSV [--outdir results]
R_ENSEMBLE_DIR is lulcc-goldens' goldens/ensemble (summary_*.csv).
"""
from __future__ import annotations

import argparse
import json
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

# Values printed in the paper (Figs. 5 and 6), one unseeded run of lulcc.
PAPER = {
    "auc_glm_Forest": 0.6423, "auc_glm_Built": 0.9384, "auc_glm_Other": 0.7276,
    "auc_builtgain": 0.6448,
}


def r_wide(d: Path) -> pd.DataFrame:
    auc = pd.read_csv(d / "summary_auc.csv")
    auc = auc[auc.model == "glm"].assign(metric=lambda x: "auc_glm_" + x.label)
    gain = pd.read_csv(d / "summary_auc_builtgain.csv").assign(metric="auc_builtgain")
    fom = pd.read_csv(d / "summary_fom_forest_built.csv")
    fom = fom.assign(metric=lambda x: "fom_forest_built_" + x.model + "_" + x.factor.astype(str))
    long = pd.concat([auc[["seed", "metric", "auc"]].rename(columns={"auc": "value"}),
                      gain[["seed", "metric", "auc"]].rename(columns={"auc": "value"}),
                      fom[["seed", "metric", "fom"]].rename(columns={"fom": "value"})])
    return long.pivot(index="seed", columns="metric", values="value")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("r_dir", type=Path)
    ap.add_argument("py_csv", type=Path)
    ap.add_argument("--outdir", type=Path, default=Path("results"))
    ap.add_argument("--acceptance", type=Path, default=Path("acceptance.toml"))
    args = ap.parse_args()
    acc = tomllib.loads(args.acceptance.read_text())["ensemble"]

    R = r_wide(args.r_dir)
    P = pd.read_csv(args.py_csv).set_index("seed")
    metrics = [m for m in R.columns if m in P.columns]

    rows = []
    for m in metrics:
        r, p = R[m].dropna().to_numpy(), P[m].dropna().to_numpy()
        se = np.sqrt(r.var(ddof=1) / len(r) + p.var(ddof=1) / len(p))
        tol = max(3 * se, acc["abs_floor"])
        diff = p.mean() - r.mean()
        ks = ks_2samp(r, p)
        row = {
            "metric": m, "n_R": len(r), "n_py": len(p),
            "mean_R": r.mean(), "sd_R": r.std(ddof=1), "min_R": r.min(), "max_R": r.max(),
            "mean_py": p.mean(), "sd_py": p.std(ddof=1),
            "diff": diff, "tolerance": tol, "means_agree": bool(abs(diff) <= tol),
            "ks_p": ks.pvalue, "distributions_agree": bool(ks.pvalue >= acc["ks_alpha"]),
        }
        row["pass"] = row["means_agree"] and row["distributions_agree"]
        if m in PAPER:
            row["paper"] = PAPER[m]
            row["paper_percentile_in_R"] = float((r <= PAPER[m]).mean() * 100)
        rows.append(row)
    res = pd.DataFrame(rows)
    args.outdir.mkdir(parents=True, exist_ok=True)
    res.to_csv(args.outdir / "ensemble_comparison.csv", index=False)
    summary = {"metrics": len(res), "passed": int(res["pass"].sum()),
               "failed": res.loc[~res["pass"], "metric"].tolist(),
               "accepted": bool(res["pass"].all())}
    (args.outdir / "ensemble_summary.json").write_text(json.dumps(summary, indent=2))

    pd.set_option("display.width", 200)
    cols = ["metric", "mean_R", "sd_R", "mean_py", "sd_py", "diff", "tolerance", "ks_p", "pass"]
    print(res[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print(json.dumps(summary, indent=2))
    _figure(R, P, args.outdir / "fig_fom_forest_built.png")
    return 0


def _figure(R: pd.DataFrame, P: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    factors = [2 ** k for k in range(1, 9)]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
    for ax, model, title in zip(axes, ("clues", "ordered"), ("CLUE-S", "Ordered")):
        cols = [f"fom_forest_built_{model}_{f}" for f in factors]
        for df, color, lab in ((R, "#1b6ca8", "lulcc (R)"), (P, "#d1495b", "DisSLUCC (Python)")):
            v = df[cols].to_numpy()
            ax.fill_between(factors, v.min(0), v.max(0), color=color, alpha=0.25, lw=0)
            ax.plot(factors, v.mean(0), "o-", color=color, ms=3, label=f"{lab}, mean and range")
        ax.set_xscale("log", base=2)
        ax.set_xticks(factors, [str(f) for f in factors])
        ax.set_title(title)
        ax.set_xlabel("Resolution (multiple of native pixel size)")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Figure of merit, forest → built")
    axes[0].legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=200)


if __name__ == "__main__":
    raise SystemExit(main())

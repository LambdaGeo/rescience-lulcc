# SPDX-License-Identifier: GPL-3.0-or-later
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
    if any(c.startswith("ab_") for c in P.columns):
        _agreement(args.r_dir, P, args.outdir)
    return 0


AB = {"a": "change simulated as persistence", "b": "change simulated correctly",
      "c": "change simulated incorrectly", "d": "persistence simulated as change"}
FACTORS = [2 ** k for k in range(1, 9)]


def r_agreement(d: Path) -> pd.DataFrame:
    """Agreement budget of forest -> built of every seed of lulcc, one row per seed."""
    rows = []
    for sd in sorted(d.glob("seed_*")):
        row = {"seed": int(sd.name.split("_")[1])}
        for model in ("clues", "ordered"):
            ab = pd.read_csv(sd / f"{model}_agreement.csv")
            ab = ab[(ab.level == "transition") & (ab.name == "Forest-Built")]
            for _, r in ab.iterrows():
                for k in AB:
                    row[f"ab_forest_built_{model}_{int(r.factor)}_{k}"] = r[k]
        rows.append(row)
    return pd.DataFrame(rows).set_index("seed")


def _agreement(r_dir: Path, P: pd.DataFrame, outdir: Path) -> None:
    """Descriptive comparison of the agreement budget (Fig. 7 of the paper).

    Not part of the acceptance criteria, which were fixed before the runs for the 20
    metrics above; the same statistics are reported for information only.
    """
    R = r_agreement(r_dir)
    rows = []
    for m in [c for c in R.columns if c in P.columns]:
        r, p = R[m].to_numpy(), P[m].to_numpy()
        se = np.sqrt(r.var(ddof=1) / len(r) + p.var(ddof=1) / len(p))
        rows.append({"metric": m, "mean_R": r.mean(), "sd_R": r.std(ddof=1),
                     "mean_py": p.mean(), "sd_py": p.std(ddof=1), "diff": p.mean() - r.mean(),
                     "diff_in_se": (p.mean() - r.mean()) / se if se > 0 else 0.0,
                     "ks_p": ks_2samp(r, p).pvalue})
    res = pd.DataFrame(rows)
    res.to_csv(outdir / "agreement_comparison.csv", index=False)
    print(f"agreement budget (descriptive): {len(res)} values, "
          f"max |diff| {res['diff'].abs().max():.2e} of the study area, "
          f"max |diff|/SE {res['diff_in_se'].abs().max():.2f}, min KS p {res['ks_p'].min():.3f}")
    _agreement_figure(R, P, outdir / "fig_agreement_forest_built.png")


def _agreement_figure(R: pd.DataFrame, P: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    colors = {"a": "#d1495b", "b": "#66a182", "c": "#edae49", "d": "#00798c"}
    x = np.arange(len(FACTORS))
    w = 0.38
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), sharey=True)
    for ax, model, title in zip(axes, ("clues", "ordered"), ("CLUE-S", "Ordered")):
        for off, df, hatch in ((-w / 2 - 0.01, R, None), (w / 2 + 0.01, P, "////")):
            bottom = np.zeros(len(FACTORS))
            for k in ("b", "a", "c", "d"):
                v = 100 * df[[f"ab_forest_built_{model}_{f}_{k}" for f in FACTORS]].mean().to_numpy()
                ax.bar(x + off, v, w, bottom=bottom, color=colors[k], hatch=hatch,
                       edgecolor="white", lw=0.4)
                bottom += v
        ax.set_xticks(x, [str(f) for f in FACTORS])
        ax.set_title(title)
        ax.set_xlabel("Resolution (multiple of native pixel size)")
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Percentage of study area")
    handles = [Patch(color=colors[k], label=AB[k]) for k in ("b", "a", "c", "d")]
    handles += [Patch(facecolor="white", edgecolor="0.3", label="left bar: lulcc (R)"),
                Patch(facecolor="white", edgecolor="0.3", hatch="////", label="right bar: replication")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.savefig(path, dpi=200)


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

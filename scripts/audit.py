"""Audit of the deterministic comparison: is the comparison what it claims to be?

Usage: python scripts/audit.py GOLDENS_DIR [--out results/audit.json] [--power]

Checks (each one fails the run if it does not hold):

  provenance   the CLUE-S that runs is DisSLUCC's published AllocationDClueSLike,
               imported from the installed package, called once per step
  inputs       from the golden only the probabilities and the demand are read; the
               simulated maps of lulcc are read only to compare
  change       lulcc's maps really change from step to step, so there is something
               to reproduce, and ours match at the same step and not at the
               neighbouring steps
  zero-jitter  the subclass that adds the perturbation, with perturbation 0, gives
               the same maps as the published component
  disslucc     the whole DisSLUCC path, with its own PotentialDLogisticRegression
               (float32) from lulcc's coefficients, also gives lulcc's maps
  controls     one change at a time (parameters, inputs) must make the maps differ:
               the comparison can fail

--power adds the ensemble negative control (about 15 minutes): 10 runs of the
replication with elasticity 0.25 must fail the criteria of acceptance.toml, and 10
runs without perturbation are reported (the criteria do not detect its absence).

Note on DisSModel: a Model attaches to the most recent Environment. Every scenario
below creates its Environment before any model; a model created while an older,
finished Environment is the current one silently never runs.
"""
from __future__ import annotations

import argparse
import inspect
import json
import sys
import tomllib
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
warnings.filterwarnings("ignore")

import disslucc  # noqa: E402
from disslucc.components.allocation import AllocationDClueSLike  # noqa: E402
from disslucc.components.demand import DemandPreComputedValues  # noqa: E402
from disslucc.components.potential import PotentialDLogisticRegression  # noqa: E402
from disslucc.schemas import LogisticRegressionSpec  # noqa: E402
from dissmodel.core import Environment  # noqa: E402

from re_lulcc import clues, data  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("goldens", type=Path)
    ap.add_argument("--out", type=Path, default=Path("results/audit.json"))
    ap.add_argument("--power", action="store_true")
    args = ap.parse_args()
    D = args.goldens / "deterministic"

    pie = data.load_pie(args.goldens / "data")
    mask = pie.mask
    gp = pd.read_csv(D / "probabilities.csv.gz")
    gd = pd.read_csv(D / "demand.csv")[data.LABELS].to_numpy()
    base_prob = {lab: gp[lab].to_numpy() for lab in data.LABELS}
    order = gp["cell"].to_numpy() - 1
    with rasterio.open(D / "clues_maps.tif") as r:
        gold = r.read().astype(np.int16)
    gold[:, ~mask] = 0
    report: dict = {}
    failures: list[str] = []

    def check(name: str, ok: bool) -> None:
        report.setdefault("checks", {})[name] = bool(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name}")
        if not ok:
            failures.append(name)

    def full(prob: dict) -> dict:
        out = {}
        for lab in data.LABELS:
            a = np.zeros(pie.shape)
            a.ravel()[order] = prob[lab]
            out[lab] = a
        return out

    def differing(maps: list[np.ndarray]) -> list[int]:
        return [int(((m != g) & mask).sum()) for m, g in zip(maps, gold)]

    def run(prob=base_prob, dmd=gd, **kw) -> list[np.ndarray]:
        b = data.make_backend(pie, pie.lu[0])
        maps = [pie.lu[0].copy()]
        clues.run_clues(b, full(prob), dmd, record=lambda t, be: maps.append(data.to_codes(be, mask)), **kw)
        return maps

    # provenance
    src = Path(inspect.getsourcefile(AllocationDClueSLike)).resolve()
    pkg = Path(disslucc.__file__).resolve().parent
    report["disslucc"] = {"version": getattr(disslucc, "__version__", None), "file": str(src)}
    check("provenance: AllocationDClueSLike comes from the installed disslucc package",
          src.is_relative_to(pkg) and not src.is_relative_to(ROOT))
    calls = {"n": 0}
    orig = AllocationDClueSLike.execute

    def counting(self):
        calls["n"] += 1
        return orig(self)

    AllocationDClueSLike.execute = counting
    try:
        ours = run()
    finally:
        AllocationDClueSLike.execute = orig
    report["published_execute_calls"] = calls["n"]
    check("provenance: published execute() called once per step (14)", calls["n"] == 14)

    # change and alignment
    changes = [int(((gold[t] != gold[t - 1]) & mask).sum()) for t in range(1, 15)]
    report["lulcc_cells_changed_per_step"] = changes
    check("change: lulcc's map changes at every step", min(changes) > 0)
    same = differing(ours)
    lag = [int(((ours[t] != gold[t - 1]) & mask).sum()) for t in range(1, 15)]
    report["differing_same_step"] = same
    check("alignment: identical at the same step", max(same) == 0)
    check("alignment: different from the previous step", min(lag) > 0)

    # zero-jitter subclass equals the published component
    b = data.make_backend(pie, pie.lu[0])
    maps = [pie.lu[0].copy()]
    env = Environment(start_time=1, end_time=14)
    dem = DemandPreComputedValues(annual_demand=[list(map(float, r)) for r in gd], land_use_types=data.LABELS)
    clues.ProbabilityPotential(backend=b, probability=full(base_prob), elasticity={l: 0.2 for l in data.LABELS})
    clues.JitteredDClueS(backend=b, demand=dem, land_use_types=data.LABELS, transition_matrix=[[[1, 1, 1]] * 3],
                         cell_area=1.0, max_difference=50, max_iteration=1000, factor_iteration=1e-6,
                         jitter=0.0, rng=np.random.default_rng(1))
    clues._Recorder(backend=b, callback=lambda t, be: maps.append(data.to_codes(be, mask)))
    env.run()
    check("zero-jitter: JitteredDClueS(jitter=0) reproduces lulcc", max(differing(maps)) == 0)

    # the whole DisSLUCC path
    co = pd.read_csv(D / "glm_coefficients.csv")
    specs = []
    for lab in data.LABELS:
        r = co[co.label == lab].set_index("term")["estimate"]
        specs.append(LogisticRegressionSpec(const=float(r["(Intercept)"]), elasticity=0.2,
                                            betas={k: float(v) for k, v in r.items() if k != "(Intercept)"}))
    b = data.make_backend(pie, pie.lu[0])
    maps = [pie.lu[0].copy()]
    env = Environment(start_time=1, end_time=14)
    dem = DemandPreComputedValues(annual_demand=[list(map(float, r)) for r in gd], land_use_types=data.LABELS)
    PotentialDLogisticRegression(backend=b, potential_data=[specs], land_use_types=data.LABELS)
    AllocationDClueSLike(backend=b, demand=dem, land_use_types=data.LABELS, transition_matrix=[[[1, 1, 1]] * 3],
                         cell_area=1.0, max_difference=50, max_iteration=1000, factor_iteration=1e-6)
    clues._Recorder(backend=b, callback=lambda t, be: maps.append(data.to_codes(be, mask)))
    env.run()
    check("disslucc: PotentialDLogisticRegression + AllocationDClueSLike reproduce lulcc",
          max(differing(maps)) == 0)

    # negative controls
    rng = np.random.default_rng(0)
    noisy = {k: v + rng.normal(0, 1e-4, v.shape) for k, v in base_prob.items()}
    shifted = gd.copy()
    shifted[1:, 1] += 100
    shifted[1:, 0] -= 100
    controls = {
        "elasticity 0.2 -> 0.25": dict(elasticity=0.25),
        "scale factor 1e-6 -> 2e-6": dict(scale=2e-6),
        "max difference 50 -> 10": dict(max_difference=10),
        "probabilities + N(0, 1e-4)": dict(prob=noisy),
        "demand +100 built, -100 forest": dict(dmd=shifted),
    }
    report["negative_controls_final_cells_differing"] = {}
    for name, kw in controls.items():
        d = differing(run(**kw))[-1]
        report["negative_controls_final_cells_differing"][name] = d
        print(f"      control {name:34s} {d:>6} cells differ")
        check(f"control detected: {name}", d > 0)

    if args.power:
        report["ensemble_power"] = ensemble_power(args.goldens, pie, check)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2))
    print(f"\n{len(failures)} failure(s)" if failures else "\nall checks passed")
    return 1 if failures else 0


def ensemble_power(goldens: Path, pie, check) -> dict:
    from scipy.stats import ks_2samp

    from compare_ensemble import r_wide
    from re_lulcc import ensemble

    R = r_wide(goldens / "ensemble")
    acc = tomllib.loads((ROOT / "acceptance.toml").read_text())["ensemble"]
    metrics = [c for c in R.columns if c.startswith("fom_forest_built_clues")]
    base = dict(ensemble.CLUES)
    out = {}
    for name, cfg in (("elasticity 0.25", dict(base, elasticity=0.25)), ("no perturbation", dict(base, jitter=0.0))):
        ensemble.CLUES = cfg
        P = pd.DataFrame([ensemble.run_seed(pie, s) for s in range(1, 11)]).set_index("seed")
        failed = []
        for m in metrics:
            r, p = R[m].to_numpy(), P[m].to_numpy()
            se = np.sqrt(r.var(ddof=1) / len(r) + p.var(ddof=1) / len(p))
            ok = abs(p.mean() - r.mean()) <= max(3 * se, acc["abs_floor"]) and ks_2samp(r, p).pvalue >= acc["ks_alpha"]
            if not ok:
                failed.append(m)
        out[name] = {"clues_fom_metrics_failing": len(failed), "of": len(metrics)}
        print(f"      ensemble variant {name:16s} fails {len(failed)}/{len(metrics)} CLUE-S metrics")
    ensemble.CLUES = base
    check("ensemble power: elasticity 0.25 fails the criteria", out["elasticity 0.25"]["clues_fom_metrics_failing"] > 0)
    return out


if __name__ == "__main__":
    raise SystemExit(main())

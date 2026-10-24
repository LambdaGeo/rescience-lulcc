# rescience-lulcc (provisional)

Replication, on DisSModel and DisSLUCC, of

> Moulds, S., Buytaert, W., & Mijic, A. (2015). An open and extensible
> framework for spatially explicit land use change modelling: the lulcc R
> package. *Geoscientific Model Development*, 8, 3215–3229.
> https://doi.org/10.5194/gmd-8-3215-2015

The target is a ReScience C submission. This repository holds the Python
replication only; the reference outputs come from
[lulcc-goldens](https://github.com/LambdaGeo/lulcc-goldens), generated with
the original R package.

**Status:** provisional. DisSLUCC 0.5.0 is used as published; the pieces it
lacks live in `src/re_lulcc/` and are meant to move into DisSLUCC later.

## What is where

| Module | Paper | Source |
|---|---|---|
| `data.py` | Sect. 3.1–3.2, data and transition tables | new |
| `suitability.py` | Sect. 3.3, binary logistic regression (statsmodels) | new |
| `demand.py` | Sect. 3.4, linear interpolation and rounding | new; rounding resolved against lulcc's `roundSum` |
| `clues.py` | Sect. 3.5.2, CLUE-S | **allocation: DisSLUCC 0.5.0 `AllocationDClueSLike`, unchanged**; provisional: `ProbabilityPotential` (probability + elasticity) and `JitteredDClueS` (jitter, warn at the iteration limit) |
| `ordered.py` | Sect. 3.5.3, ordered allocation (Fuchs et al., 2013) | provisional, see below |
| `threemap.py` | Sect. 3.6, three-map comparison, agreement budget, figure of merit (Pontius et al., 2011) | provisional, see below |

DisSLUCC's CLUE-S came from LuccME (INPE), not from lulcc. The two
implementations were written independently from Verburg et al. (2002).

Points the paper leaves open, resolved by reading lulcc's source (to be
declared in the article):
- **Ordered allocation:** which cells are candidates, the temporary class
  (-1) for released cells, the cap of 1000 random draws, ties kept in cell
  order.
- **Three-map comparison:** the ((n+1)², n+1) table layout, aggregation
  from the top-left corner with partial edge blocks weighted by their valid
  cells, and empty denominators (0/0) dropped from the sums.
- **Demand:** rows rescaled to the number of cells, floored, and the
  remaining cells given to the largest remainders (later column on ties).

## Result against the deterministic golden

`make check GOLDENS=...` (24 Oct 2026, goldens from the smoke-test
environment of lulcc-goldens, not yet the Docker release):

| Stage | Result |
|---|---|
| Observed transition tables | identical |
| glm coefficients (same training cells) | max. relative difference 1.3 × 10⁻⁹ |
| Probabilities from R's coefficients | max. difference 9 × 10⁻¹⁴ |
| Probabilities from our coefficients | max. difference 2.8 × 10⁻¹⁰ |
| Demand, t = 0..14 | identical |
| **CLUE-S (DisSLUCC 0.5.0, jitter 0), 14 steps** | **0 cells differ at every step** |
| Ordered allocation (deterministic), 14 steps | 0 cells differ at every step |
| Three-map tables, 8 resolutions | max. difference 1.1 × 10⁻¹⁶ |
| Figure of merit, forest → built | max. difference 2.2 × 10⁻¹⁶ |
| Whole chain from our own glm | final map identical |

The CLUE-S result also holds with the probabilities rounded to float32.

## Still to do

1. Ensemble: our own partition and jitter over 30 seeds against the 30
   seeds of the golden (AUCs, Figs. 5–8), with criteria written in
   `acceptance.toml` before the comparison.
2. ROC/AUC (sklearn) for Figs. 5–6, and rpart and random forest models.
3. Figures side by side with the paper's.
4. `requirements.lock`, registered run, ReScience article.
5. Move `JitteredDClueS`, `ordered.py` and `threemap.py` into DisSLUCC.

## Use

```bash
make env
make test
make check GOLDENS=../lulcc-goldens/goldens
```

Licence: MIT.

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

## Result of the ensembles

Criteria in `acceptance.toml`, committed (932d9ed) before either ensemble
ran. 30 seeds of lulcc (R) against 30 seeds of the replication, 20 metrics
(glm AUCs of Figs. 5–6; figure of merit forest → built of Fig. 8 at 8
resolutions, for CLUE-S and ordered allocation): **20 of 20 pass** (means
within 3 standard errors and Kolmogorov–Smirnov p ≥ 0.01). Details in
`results/ensemble_comparison.csv`; figure in
`results/fig_fom_forest_built.png`.

Where the paper's printed AUCs (one unseeded run) fall in lulcc's own 30
runs: built 0.9384 (70th percentile), gain of built 0.6448 (63rd), other
0.7276 (7th), forest 0.6423 (just below the minimum of the 30 runs,
0.6425). The paper used lulcc 1.0 in 2015; the goldens use 1.0.4 on
current R.

## Still to do

1. ~~Regenerate the goldens in the Docker environment and rerun `make check`
   and `make compare` against them~~ (done 25 Oct 2026: same results).
2. rpart and random forest models of Fig. 5 (descriptive only).
3. Figures side by side with the paper's.
4. `requirements.lock`, registered run, ReScience article.
5. Move `JitteredDClueS`, `ordered.py` and `threemap.py` into DisSLUCC.

## Audit

`make audit` (about 3 minutes) checks that the comparison is what it claims: the
CLUE-S that runs is DisSLUCC's published component, called once per step; only
suitability and demand are read from the reference; lulcc's maps change at every step
and ours match at the same step only; the perturbation subclass with perturbation 0
equals the published component; DisSLUCC's own logistic potential also reproduces
lulcc; and five negative controls (elasticity, scale factor, tolerance, noise on the
suitability, demand) make the maps differ. `make audit-power` adds the ensemble
negative control (about 15 minutes more). Results in `results/audit.json`.

## Use

```bash
make env
make test
make check GOLDENS=../lulcc-goldens/goldens
```

## Licence

Copyright (C) 2026 Sérgio Souza Costa, Denilson da Silva Bezerra and José Magno Pinheiro
Alves (UFMA).

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. It is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; see `LICENSE`.

GPL-3.0-or-later, the licence family of lulcc (GPL), because parts of this code (the
ordered allocation, the three-map comparison, the demand rounding and the training
sample) follow lulcc's source code step by step. Dependencies keep their own licences
(DisSModel and DisSLUCC 0.5.0: MIT). The article itself is published under CC BY 4.0
by ReScience C.

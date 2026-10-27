# rescience-lulcc

Replication, on DisSModel and DisSLUCC, of

> Moulds, S., Buytaert, W., & Mijic, A. (2015). An open and extensible
> framework for spatially explicit land use change modelling: the lulcc R
> package. *Geoscientific Model Development*, 8, 3215–3229.
> https://doi.org/10.5194/gmd-8-3215-2015

The target is a ReScience C submission. This repository holds the Python
replication only; the reference outputs come from
[lulcc-goldens](https://github.com/LambdaGeo/lulcc-goldens), generated with
the original R package.

**Status:** v0.1.0, the version described in the article (`article/article.pdf`).
DisSLUCC 0.5.0 is used as published; the pieces it lacks live in `src/re_lulcc/` and
are planned for DisSLUCC 0.6.0 (then v0.2.0 of this repository).

## What is where

Origins as in Fig. 1 and Section 2 of the article: *DisSLUCC* (used unchanged),
*written* for the replication without reference to lulcc's code, or *port* of lulcc's
code (followed step by step where the paper leaves a choice open).

| Module | Paper | Origin |
|---|---|---|
| `data.py` | Sect. 3.1–3.2, data and transition tables | written |
| `suitability.py` | Sect. 3.3, binary logistic regression (statsmodels) | written |
| `demand.py` | Sect. 3.4, linear interpolation and rounding | interpolation written; rounding a port of `roundSum` |
| `clues.py` | Sect. 3.5.2, CLUE-S | **allocation: DisSLUCC 0.5.0 `AllocationDClueSLike`, unchanged**; `ProbabilityPotential` (probability + elasticity) written; `JitteredDClueS` (random perturbation, warn at the iteration limit) a port |
| `ordered.py` | Sect. 3.5.3, ordered allocation (Fuchs et al., 2013) | port |
| `threemap.py` | Sect. 3.6, three-map comparison, agreement budget, figure of merit (Pontius et al., 2011) | port |
| `ensemble.py` | one run of the paper's experiment per seed | AUC written; training sample a port of `caret::createDataPartition` |

DisSLUCC's CLUE-S came from LuccME (INPE), not from lulcc. The two
implementations were written independently from Verburg et al. (2002).

Points the paper leaves open, resolved by reading lulcc's source (Section 4 of the
article):
- **Ordered allocation:** which cells are candidates, the temporary class
  (-1) for released cells, the cap of 1000 random draws, ties kept in cell
  order.
- **Three-map comparison:** the ((n+1)², n+1) table layout, aggregation
  from the top-left corner with partial edge blocks weighted by their valid
  cells, and empty denominators (0/0) dropped from the sums.
- **Demand:** rows rescaled to the number of cells, floored, and the
  remaining cells given to the largest remainders (later column on ties).

## Result against the deterministic golden

`make check GOLDENS=...` against the goldens of lulcc-goldens v0.1.0 (Docker image,
R 4.5.1):

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

The agreement budget of forest → built (Fig. 7 of the paper) was added after the
ensembles had run, so it is outside the criteria and reported for information: over
its 64 values, means within 1.02 standard errors and KS p ≥ 0.24
(`results/agreement_comparison.csv`, `results/fig_agreement_forest_built.png`).

Where the paper's printed AUCs (one unseeded run) fall in lulcc's own 30
runs: built 0.9384 (70th percentile), gain of built 0.6448 (63rd), other
0.7276 (7th), forest 0.6423 (just below the minimum of the 30 runs,
0.6425). The paper used lulcc 1.0 in 2015; the goldens use 1.0.4 on
current R.

## Still to do

1. rpart and random forest models of Fig. 5 (descriptive only).
2. Move `JitteredDClueS`, `ordered.py` and `threemap.py` into DisSLUCC (DisSLUCC 0.6.0,
   then v0.2.0 of this repository).

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

With [lulcc-goldens](https://github.com/LambdaGeo/lulcc-goldens) cloned next to this
repository (otherwise pass `GOLDENS=path/to/lulcc-goldens/goldens`):

```bash
make env        # .venv with the pinned packages (requirements.lock)
make test       # unit tests, seconds
make check      # deterministic level, Table 1 (~1 min)
make audit      # checks of the comparison, Section 8.1 (~3 min)
make maps       # Fig. 2
make ensemble   # 30 runs of the replication (~15 min)
make compare    # acceptance criteria, Table 2, Figs. 3 and 4
make figures    # copy the figures of results/ into article/figures/
```

## Building the article

The article is in `article/` (ReScience C template, XeLaTeX and biber). On Ubuntu or
Debian:

```bash
sudo apt install make latexmk texlive-xetex texlive-latex-extra \
    texlive-fonts-recommended texlive-bibtex-extra biber python3-yaml
sudo locale-gen en_US.UTF-8
cd article && make
```

Build it outside the replication's virtual environment (`deactivate` first): the
template's `yaml-to-latex.py` needs PyYAML, which the system package `python3-yaml`
provides and `.venv` does not. `make` also rebuilds Fig. 1
(`article/figures/fig_method.tex`, TikZ) when its source changes; Figs. 2–4 come from
`results/` through `make figures`. After installing missing packages, run `make clean`
before `make`.

## Licence

Copyright (C) 2026 Sérgio Souza Costa, Denilson da Silva Bezerra, José Magno Pinheiro Alves,
Hissa Barbara Oliveira, Mateus Dutra Vale and Samylli Kalei Silva Zhao (UFMA).

This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. It is
distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; see `LICENSE`.

GPL-3.0-or-later, the licence family of lulcc (GPL), because parts of this code (the
ordered allocation, the three-map comparison, the demand rounding and the training
sample) follow lulcc's source code step by step. Dependencies keep their own licences
(DisSModel and DisSLUCC 0.5.0: MIT). The article itself is published under CC BY 4.0
by ReScience C.

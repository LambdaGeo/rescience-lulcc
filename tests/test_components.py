import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from re_lulcc import demand, ordered, threemap  # noqa: E402


def test_round_sum_rows_sum_to_ncell_and_largest_remainders_win():
    x = np.array([[1.4, 1.4, 1.2]])  # sums to 4; floors 1,1,1; one extra cell
    out = demand.round_sum(x, 4)
    assert out.sum() == 4
    # equal remainders: the later column wins, as in lulcc's roundSum (stable sort, take the last)
    assert out.tolist() == [[1, 2, 1]]


def test_interpolation_hits_observed_totals():
    totals = np.array([[10, 5, 5], [8, 7, 5], [4, 10, 6]])
    d = demand.interpolated_demand([0, 6, 14], totals, np.arange(15), 20)
    assert (d[[0, 6, 14]] == totals).all()
    assert (d.sum(axis=1) == 20).all()


def test_ordered_meets_demand_deterministically():
    rng = np.random.default_rng(0)
    prob = rng.uniform(size=(100, 3))
    lu0 = rng.integers(1, 4, size=100)
    target = np.array([30, 40, 30])
    lu1 = ordered.ordered_step(prob, lu0, target, [1, 2, 3], [2, 1, 3], stochastic=False)
    assert [(lu1 == c).sum() for c in (1, 2, 3)][1] == 40


def test_identical_maps_give_perfect_figure_of_merit():
    rng = np.random.default_rng(1)
    x = rng.integers(1, 4, size=(20, 20))
    x1 = x.copy()
    x1[:5] = 2
    tabs = threemap.tables(x, x1, x1, [1, 2], [1, 2, 3])
    for tab in tabs:
        fom = threemap.figure_of_merit(tab, 3)
        assert np.isclose(fom["overall"], 1.0)

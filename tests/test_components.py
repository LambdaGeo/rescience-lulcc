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


def test_jitter_subclass_at_zero_equals_published_component():
    """JitteredDClueS replaces the iteration loop of DisSLUCC's AllocationDClueSLike to add
    the perturbation; with perturbation 0 it must give exactly the published result."""
    from disslucc.components.allocation import AllocationDClueSLike
    from disslucc.components.demand import DemandPreComputedValues
    from dissmodel.core import Environment
    from dissmodel.geo import RasterBackend

    from re_lulcc import clues
    from re_lulcc.data import LABELS

    rng = np.random.default_rng(3)
    shape = (40, 50)
    lu0 = rng.integers(0, 3, size=shape)
    prob = {lab: rng.uniform(size=shape) for lab in LABELS}
    totals = np.array([(lu0 == k).sum() for k in range(3)])
    demand = [totals.tolist(), (totals + [-60, 80, -20]).tolist(), (totals + [-120, 150, -30]).tolist()]

    def run(cls, **kw):
        b = RasterBackend(shape=shape)
        for k, lab in enumerate(LABELS):
            b.set(lab, (lu0 == k).astype(np.float32))
        env = Environment(start_time=1, end_time=2)
        dem = DemandPreComputedValues(annual_demand=[list(map(float, r)) for r in demand], land_use_types=LABELS)
        clues.ProbabilityPotential(backend=b, probability=prob, elasticity={lab: 0.2 for lab in LABELS})
        cls(backend=b, demand=dem, land_use_types=LABELS, transition_matrix=[[[1, 1, 1]] * 3], cell_area=1.0,
            max_difference=10, max_iteration=5000, factor_iteration=1e-4, **kw)
        env.run()
        return np.stack([b.get(lab) for lab in LABELS])

    published = run(AllocationDClueSLike)
    ours = run(clues.JitteredDClueS, jitter=0.0, rng=np.random.default_rng(1))
    assert np.array_equal(published, ours)
    assert not np.array_equal(published[:, :, :], np.stack([(lu0 == k) for k in range(3)]).astype(np.float32))

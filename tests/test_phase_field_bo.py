import logging

import numpy as np
import pandas as pd
import pytest
from gpytorch.kernels import MaternKernel, RBFKernel

from phasebo.phase_field_bo import PhaseFieldBO
from phasebo.utils.other import set_seeds

NEW = ['Li9 S2 Cl5', 'Li6 S1 Cl4', 'Li11 S3 Cl5']
COMPUTED = ['Li3 S1 Cl1', 'Li4 S1 Cl2']


def make_bo(li_s_cl, **kwargs):
    compositions, references, ions = li_s_cl
    set_seeds(0)
    return PhaseFieldBO(compositions, references, ions, exclude_zeros=True, logger=logging.getLogger('test_logger'),
                        **kwargs)


def rows(X):
    return [tuple(x) for x in X]


def names(bo, X):
    return {bo.next_list[bo.fcsym(x)] for x in X}


def test_unevaluated_drops_duplicates_and_evaluated_points():
    choices = np.array([[0.1, 0.2], [0.3, 0.4], [0.1, 0.2], [0.5, 0.6]])
    assert rows(PhaseFieldBO.unevaluated(choices, np.array([[0.3, 0.4]]))) == [(0.1, 0.2), (0.5, 0.6)]


@pytest.mark.parametrize('acquisition', ['qlogei', 'ts'])
def test_suggest_picks_distinct_new_candidates(li_s_cl, acquisition):
    bo = make_bo(li_s_cl, mode='suggest', next_formulas=NEW + COMPUTED, batch=2, acquisition=acquisition)
    assert len(set(rows(bo.next))) == 2
    assert names(bo, bo.next) <= set(NEW)


@pytest.mark.parametrize('acquisition', ['qlogei', 'ts'])
def test_suggest_batch_shrinks_to_the_candidates_left(li_s_cl, acquisition):
    bo = make_bo(li_s_cl, mode='suggest', next_formulas=NEW + COMPUTED, batch=5, acquisition=acquisition)
    assert names(bo, bo.next) == set(NEW)


@pytest.mark.parametrize('mode, options', [('suggest', {'next_formulas': NEW}), ('path', {'n_seeds': 2, 'max_iter': 1})])
def test_references_are_not_training_data(li_s_cl, mode, options):
    bo = make_bo(li_s_cl, mode=mode, batch=2, **options)
    references = {tuple(bo.dic[name][1]) for name in bo.references}
    assert not references & set(rows(bo.model.train_inputs[0].numpy()))


def test_suggest_skips_the_references(li_s_cl):
    bo = make_bo(li_s_cl, mode='suggest', next_formulas=NEW + ['Li2 S1', 'Li4 Cl4'], batch=5)
    assert names(bo, bo.next) == set(NEW)


def test_suggest_batch_is_empty_when_every_candidate_is_computed(li_s_cl):
    bo = make_bo(li_s_cl, mode='suggest', next_formulas=COMPUTED, batch=2)
    assert len(bo.next) == 0


@pytest.mark.parametrize('acquisition', ['qlogei', 'ts'])
def test_path_evaluates_new_points_at_their_lowest_energy(li_s_cl, acquisition):
    bo = make_bo(li_s_cl, mode='path', n_seeds=2, max_iter=2, batch=2, acquisition=acquisition)
    seeds, evals = rows(bo.X[:2]), rows(bo.X[2:])
    assert len(set(evals)) == len(evals) == 4
    assert not set(evals) & set(seeds)
    assert np.allclose(bo.Y[2:, 0], [bo.f(x) for x in bo.X[2:]])


def test_path_stops_when_every_candidate_is_evaluated(li_s_cl):
    bo = make_bo(li_s_cl, mode='path', n_seeds=2, max_iter=10, batch=2)
    # candidates_fc excludes the seeds; a seed can share a point with a candidate
    assert set(rows(bo.X[2:])) == set(rows(bo.candidates_fc)) - set(rows(bo.X[:2]))


def test_uncertainty_csv(li_s_cl, tmp_path):
    bo = make_bo(li_s_cl, mode='suggest', next_formulas=NEW, batch=2, output_dir=tmp_path)
    bo.get_uncertainty()
    df = pd.read_csv(tmp_path / 'posterior.csv')
    assert list(df.columns) == ['Candidates', 'Posterior mean (meV/atom)', 'Posterior std (meV/atom)']
    assert sorted(df['Candidates']) == sorted(NEW)
    assert (df['Posterior std (meV/atom)'] > 0).all()


def test_path_csv(li_s_cl, tmp_path):
    bo = make_bo(li_s_cl, mode='path', n_seeds=2, max_iter=2, batch=2, output_dir=tmp_path)
    bo.print_results()
    df = pd.read_csv(tmp_path / 'bo_path.csv')
    assert list(df['composition'][:2]) == list(bo.seeds)
    assert list(df['seed']) == [True] * 2 + [False] * 4
    assert list(df['energy (meV/atom)']) == list(bo.Y.ravel().round(2))


@pytest.mark.parametrize('kernel, cls', [('matern', MaternKernel), ('rbf', RBFKernel)])
def test_kernel_option(li_s_cl, kernel, cls):
    bo = make_bo(li_s_cl, mode='suggest', next_formulas=NEW, batch=2, kernel=kernel)
    assert [type(k) for k in bo.model.covar_module.modules() if isinstance(k, (MaternKernel, RBFKernel))] == [cls]


def test_unsupported_kernel(li_s_cl):
    with pytest.raises(ValueError, match='Unsupported kernel'):
        make_bo(li_s_cl, mode='suggest', next_formulas=NEW, kernel='linear')


def test_unsupported_acquisition(li_s_cl):
    with pytest.raises(ValueError, match='Unsupported acquisition'):
        make_bo(li_s_cl, mode='suggest', next_formulas=NEW, acquisition='ucb')

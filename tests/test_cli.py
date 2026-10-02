import os
import subprocess
import sys

import pandas as pd
import pytest
from hydra import compose, initialize_config_dir
from hydra.errors import ConfigCompositionException
from omegaconf import OmegaConf
from pydantic import ValidationError

from phasebo.__main__ import CONFIG_DIR
from phasebo.schema import PhaseBOConfig

NEW = ['Li6 Sn1 S5', 'Li3 Sn1 S2 Cl3', 'Li1 Sn1 S1 Cl3']
PATH = ['mode=path', 'mode.n_seeds=2', 'mode.max_iter=2']


def compose_config(overrides=()):
    with initialize_config_dir(config_dir=str(CONFIG_DIR), version_base='1.3'):
        return compose(config_name='config', overrides=list(overrides))


def validate(tmp_path, overrides=()):
    """The default config, with an input file that exists, plus overrides."""
    (tmp_path / 'field.csv').touch()
    return PhaseBOConfig.from_hydra(compose_config([f'system.inputfile={tmp_path / "field.csv"}', *overrides]))


def test_default_config_is_valid_from_any_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    PhaseBOConfig.from_hydra(compose_config())


@pytest.mark.parametrize('mode', ['suggest', 'path', 'generate'])
def test_each_mode_is_valid(tmp_path, mode):
    assert validate(tmp_path, [f'mode={mode}']).mode.name == mode


@pytest.mark.parametrize('override', [
    'system.exludefile=x',  # misspelt key
    'mode=foo',
    'mode.n_seeds=10',      # path-only option
])
def test_config_rejects(override):
    with pytest.raises(ConfigCompositionException):
        compose_config([override])


@pytest.mark.parametrize('overrides', [
    ['bo.acquisition=ucb'],
    ['bo.kernel=linear'],
    ['mode=path', 'mode.n_seeds=0'],
    ['mode=path', 'mode.max_iter=abc'],
    ['system.inputfile=missing.csv'],
    ['paths.root_dir=missing'],
    ['~system.limits.Cl'],          # limits without an ion
    ['+system.exludefile=x'],       # misspelt key in a system file
])
def test_validation_rejects(tmp_path, overrides):
    with pytest.raises(ValidationError):
        validate(tmp_path, overrides)


@pytest.fixture
def config_dir(li_sn_s_cl, tmp_path):
    """A config dir with system=synthetic: the li_sn_s_cl field, NEW as its candidates."""
    compositions, references, ions = li_sn_s_cl
    pd.DataFrame(compositions).to_csv(tmp_path / 'field.csv', index=False)
    pd.DataFrame({'formula': NEW}).to_csv(tmp_path / 'candidates.csv', index=False)
    (tmp_path / 'config' / 'system').mkdir(parents=True)
    OmegaConf.save({'inputfile': str(tmp_path / 'field.csv'), 'compositionfile': str(tmp_path / 'candidates.csv'),
                    'excludefile': None, 'reference_index': len(compositions) - len(references), 'ions': ions,
                    'limits': {el: [0, 7] for el in ions}, 'N_atom': 8},
                   tmp_path / 'config' / 'system' / 'synthetic.yaml')
    return tmp_path / 'config'


def phasebo(tmp_path, config_dir, *overrides, multirun=False, returncode=0):
    """Run from an empty directory, writing to tmp_path/outputs or tmp_path/multirun."""
    flags = ['-cd', str(config_dir)] + (['-m'] if multirun else [])
    (tmp_path / 'cwd').mkdir(exist_ok=True)
    result = subprocess.run(
        [sys.executable, '-m', 'phasebo', *flags, 'system=synthetic', 'bo.batch_size=2',
         f'paths.output_dir={tmp_path}/outputs', f'paths.multirun_dir={tmp_path}/multirun', *overrides],
        cwd=tmp_path / 'cwd', env={**os.environ, 'MPLBACKEND': 'Agg'}, capture_output=True, text=True)
    assert result.returncode == returncode, result.stderr
    return result


@pytest.mark.parametrize('mode, outputs', [
    ('generate', ['candidates_list.csv']),
    ('suggest', ['posterior.csv']),
    ('path', ['BO_Path_in_Li-Sn-S-Cl.txt', 'convergence.png']),
])
def test_run_writes_outputs_to_its_directory(tmp_path, config_dir, mode, outputs):
    phasebo(tmp_path, config_dir, *(PATH if mode == 'path' else [f'mode={mode}']), 'show_plots=false')
    (run_dir,) = tmp_path.glob(f'outputs/synthetic/*/*_{mode}')
    for name in outputs + ['convex_hull.png', 'phasebo.log', '.hydra/config.yaml']:
        assert (run_dir / name).exists(), name


def test_model_options_reach_the_model(tmp_path, config_dir):
    phasebo(tmp_path, config_dir, 'mode=suggest', 'bo.kernel=rbf', 'bo.acquisition=ts', 'show_plots=false')
    (log,) = tmp_path.glob('outputs/synthetic/*/*_suggest/phasebo.log')
    assert 'Model: rbf kernel, ts acquisition, batch 2' in log.read_text()


def test_suggest_ranks_the_compositionfile_candidates(tmp_path, config_dir):
    phasebo(tmp_path, config_dir, 'mode=suggest', 'show_plots=false')
    (posterior,) = tmp_path.glob('outputs/synthetic/*/*_suggest/posterior.csv')
    assert sorted(pd.read_csv(posterior)['Candidates']) == sorted(NEW)


def test_multirun_seeds_each_job(tmp_path, config_dir):
    phasebo(tmp_path, config_dir, *PATH, 'seed=0,0,1', multirun=True)
    paths = [p.read_text() for p in sorted(tmp_path.glob('multirun/synthetic/*/*/*/BO_Path_in_*.txt'))]
    assert len(paths) == 3
    assert paths[0] == paths[1] != paths[2]


def test_experiment_summarises_each_job(tmp_path, config_dir):
    phasebo(tmp_path, config_dir, 'experiment=compare_acquisitions', 'mode.n_seeds=2', 'mode.max_iter=2', 'seed=0,1')
    (summary,) = tmp_path.glob('multirun/synthetic/*/*/summary.csv')
    df = pd.read_csv(summary)
    assert sorted(zip(df['bo.acquisition'], df['seed'])) == [('qlogei', 0), ('qlogei', 1), ('ts', 0), ('ts', 1)]
    assert (df['evaluations'] == 4).all()   # mode.max_iter x bo.batch_size


def test_invalid_config_exits_with_an_error(tmp_path, config_dir):
    result = phasebo(tmp_path, config_dir, 'bo.acquisition=ucb', returncode=1)
    assert 'Invalid configuration' in result.stdout

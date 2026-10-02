import pandas as pd
from omegaconf import OmegaConf

from phasebo.callbacks import PathSummary


def write_job(sweep_dir, num, task, energies):
    """A job directory as a 'path' run leaves it; the first composition is the seed."""
    job = sweep_dir / f'{num}_job'
    (job / '.hydra').mkdir(parents=True)
    OmegaConf.save({'hydra': {'job': {'num': num}, 'overrides': {'task': task}}}, job / '.hydra' / 'hydra.yaml')
    pd.DataFrame({'composition': [f'C{i}' for i in range(len(energies))], 'energy (meV/atom)': energies,
                  'seed': [True] + [False] * (len(energies) - 1)}).to_csv(job / 'bo_path.csv', index=False)


def summarise(sweep_dir):
    PathSummary().on_multirun_end(OmegaConf.create({'hydra': {'sweep': {'dir': str(sweep_dir)}}}))


def test_summary_row_per_job(tmp_path):
    write_job(tmp_path, 10, ['mode=path', 'seed=1'], [0.0, 3.0, 6.0, 9.0])
    write_job(tmp_path, 2, ['mode=path', 'seed=0'], [0.0, 4.0, -0.0, 11.0])
    summarise(tmp_path)
    assert pd.read_csv(tmp_path / 'summary.csv').to_dict('list') == {
        'job': [2, 10],             # by number, not as text
        'seed': [0, 1],             # swept; mode=path is the same in both
        'evaluations': [3, 3],      # the seed composition is not an evaluation
        'stable found': [1, 0],
        'best (meV/atom)': [0.0, 3.0],
        'mean (meV/atom)': [5.0, 6.0],
    }
    assert '-0.0' not in (tmp_path / 'summary.csv').read_text()


def test_no_summary_without_path_jobs(tmp_path):
    summarise(tmp_path)
    assert not (tmp_path / 'summary.csv').exists()

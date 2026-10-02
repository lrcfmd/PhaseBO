import logging
from pathlib import Path

import pandas as pd
from hydra.experimental.callback import Callback
from omegaconf import DictConfig, OmegaConf

log = logging.getLogger(__name__)


class PathSummary(Callback):
    """After a multirun, tabulate the BO path of every 'path' job in <sweep dir>/summary.csv."""

    def on_multirun_end(self, config: DictConfig, **kwargs) -> None:
        sweep_dir = Path(config.hydra.sweep.dir)
        jobs, overrides, results = [], [], []
        for bo_path in sweep_dir.rglob('bo_path.csv'):
            job = OmegaConf.load(bo_path.parent / '.hydra' / 'hydra.yaml').hydra
            evaluated = pd.read_csv(bo_path).query('not seed')['energy (meV/atom)']
            jobs.append(job.job.num)
            overrides.append({key: value for key, _, value in (o.partition('=') for o in job.overrides.task)})
            results.append({
                'evaluations': len(evaluated),
                'stable found': int((evaluated == 0).sum()),   # bo_path.csv rounds to 0.01
                'best (meV/atom)': evaluated.min() + 0.0,      # no -0.0
                'mean (meV/atom)': round(evaluated.mean(), 2),
            })
        if not jobs:
            return

        overrides = pd.DataFrame(overrides)
        swept = overrides.loc[:, overrides.nunique(dropna=False) > 1]
        summary = pd.concat([pd.Series(jobs, name='job'), swept, pd.DataFrame(results)], axis=1).sort_values('job')
        summary.to_csv(sweep_dir / 'summary.csv', index=False)
        log.info(f"BO path summary of {len(jobs)} jobs written to {sweep_dir / 'summary.csv'}\n{summary.to_string(index=False)}")

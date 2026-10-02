import logging
import os
from typing import Dict, Optional, List

import hydra
import pandas as pd
import rootutils
from hydra.core.hydra_config import HydraConfig
from hydra.types import RunMode
from omegaconf import DictConfig, OmegaConf
from pydantic import ValidationError

from phasebo.phase_field_bo import PhaseFieldBO
from phasebo.schema import PhaseBOConfig
from phasebo.utils.other import set_seeds

# exports PROJECT_ROOT for config/paths
ROOT = rootutils.setup_root(__file__, indicator=".project-root", pythonpath=True)
CONFIG_DIR = ROOT / "config"

logger = logging.getLogger("phasebo")


def save_plot(plt, path: str, show: bool) -> None:
    plt.savefig(path, dpi=150, bbox_inches='tight')
    if show:
        plt.show()
    plt.close('all')


def run(
    compositions,
    references,
    ions: Dict[str, int],
    mode: str,
    Ntot: int,
    logger,
    seeds_type: str = 'random',
    n_seeds: int = 9,
    max_iter: int = 10,
    output_dir: str = '.',
    show_plots: bool = True,
    batch_size: int = 4,
    acquisition: str = 'qlogei',
    kernel: str = 'matern',
    disect: int = 3,
    limits: Optional[Dict[str, List[int]]] = None,
    next_formulas: Optional[List[str]] = None,
    exceptions: Optional[List[str]] = None,
    allow_negative: bool = False
) -> PhaseFieldBO:
    """Main BO run function"""
    bopt = PhaseFieldBO(
        compositions=compositions,
        references=references,
        ions=ions,
        mode=mode,
        seeds_type=seeds_type,
        n_seeds=n_seeds,
        exclude_zeros=True,
        disect=disect,
        Ntot=Ntot,
        limits=limits,
        max_iter=max_iter,
        next_formulas=next_formulas,
        batch=batch_size,
        acquisition=acquisition,
        kernel=kernel,
        exceptions=exceptions,
        allow_negative=allow_negative,
        logger=logger,
        output_dir=output_dir
    )

    save_plot(bopt.plot_convex(), os.path.join(output_dir, 'convex_hull.png'), show_plots)

    if mode == 'path':
        save_plot(bopt.plot_convergence(), os.path.join(output_dir, 'convergence.png'), show_plots)
        bopt.print_results()
    elif mode == 'suggest':
        bopt.print_results()
        bopt.get_uncertainty()

    return bopt


def read_formulas(path: Optional[str]) -> Optional[List[str]]:
    """First column of a CSV file of formulas, or None if no readable file is given."""
    if not path:
        return None
    try:
        return [i[0] for i in pd.read_csv(path).values]
    except Exception as ex:
        logger.info(f"Not using {path}: {ex}")
        return None


@hydra.main(version_base="1.3", config_path=str(CONFIG_DIR), config_name="config")
def main(cfg: DictConfig) -> None:
    try:
        config = PhaseBOConfig.from_hydra(cfg)
    except ValidationError as ex:
        logger.error(f"Invalid configuration: {ex}")
        raise SystemExit(1)
    hydra_cfg = HydraConfig.get()

    logger.info("========== CONFIGURATION ==========")
    for line in OmegaConf.to_yaml(cfg, resolve=True).splitlines():
        logger.info(line)
    logger.info("===================================")

    if config.seed is not None:
        set_seeds(config.seed)

    system = config.system
    df = pd.read_csv(system.inputfile, header=0)

    run(
        compositions=df.values,
        references=df.values[system.reference_index:],
        ions=system.ions,
        mode=config.mode.name,
        **config.mode.model_dump(exclude={'name'}),
        Ntot=system.N_atom,
        logger=logger,
        output_dir=hydra_cfg.runtime.output_dir,
        show_plots=config.show_plots and hydra_cfg.mode == RunMode.RUN,
        batch_size=config.bo.batch_size,
        acquisition=config.bo.acquisition,
        kernel=config.bo.kernel,
        limits=system.limits,
        next_formulas=read_formulas(system.compositionfile),
        exceptions=read_formulas(system.excludefile),
        allow_negative=False
    )
    logger.info(f"Outputs written to {hydra_cfg.runtime.output_dir}")


if __name__ == "__main__":
    main()

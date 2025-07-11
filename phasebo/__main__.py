import argparse
import yaml
import pandas as pd
import time
from typing import Dict, Optional, List, Any

from phasebo.phase_field_bo import PhaseFieldBO
from phasebo.logger import get_logger

def run(
    compositions,
    references,
    ions: Dict[str, int],
    mode: str,
    Ntot: int,
    seeds_type: str,
    n_seeds: int,
    max_iter: int,
    log_name: str,
    logger,
    batch_size: int = 4,
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
        Ntot=Ntot,
        limits=limits,
        max_iter=max_iter,
        next_formulas=next_formulas,
        batch=batch_size,
        exceptions=exceptions,
        allow_negative=allow_negative,
        logger=logger
    )

    convex = bopt.plot_convex()
    convex.show()

    if mode == 'path':
        bopt.bo.plot_convergence()
        bopt.bo.plot_acquisition()
        bopt.print_results()
    elif mode == 'suggest':
        bopt.print_results()
        bopt.get_uncertainty()

    return bopt

def main():
    parser = argparse.ArgumentParser(description="Run phasebo with a specified YAML configuration file.")
    parser.add_argument(
        "--config",
        type=str,
        default="input_config.yaml",
        help="Path to YAML config file (default: input_config.yaml)"
    )
    args = parser.parse_args()

    with open(args.config, "r") as f:
        cfg: Dict[str, Any] = yaml.safe_load(f)

    df = pd.read_csv(cfg["inputfile"], header=0)
    compositions = df.values
    references = df.values[cfg["reference_index"]:]

    next_formulas = None
    if "compositionfile" in cfg and cfg["compositionfile"]:
        try:
            next_formulas = [i[0] for i in pd.read_csv(cfg["compositionfile"]).values]
        except Exception:
            next_formulas = None

    exceptions = None
    if "excludefile" in cfg and cfg["excludefile"]:
        try:
            exceptions = [i[0] for i in pd.read_csv(cfg["excludefile"]).values]
        except Exception:
            exceptions = None

    # Build log file path with timestamp
    timestamp = time.strftime('%b-%d-%Y_%H%M', time.localtime())
    log_path = f"{cfg['log']}-{timestamp}.log"

    logger = get_logger("phasebo", log_file=log_path)

    # Echo config nicely
    logger.info("========== CONFIGURATION ==========")
    for k, v in cfg.items():
        logger.info(f"{k}: {v}")
    logger.info("===================================")

    run(
        compositions=compositions,
        references=references,
        ions=cfg["ions"],
        mode=cfg["mode"],
        Ntot=cfg["N_atom"],
        seeds_type=cfg["seeds_type"],
        n_seeds=cfg["n_seeds"],
        max_iter=cfg["max_iter"],
        log_name=cfg["log"],
        logger=logger,
        limits=cfg.get("limits"),
        next_formulas=next_formulas,
        exceptions=exceptions,
        allow_negative=False
    )

if __name__ == "__main__":
    main()

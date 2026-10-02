from pathlib import Path
from typing import Annotated, Dict, Literal, Optional, Tuple, Union

from omegaconf import DictConfig, OmegaConf
from pydantic import BaseModel, ConfigDict, DirectoryPath, Field, FilePath, NonNegativeInt, PositiveInt, model_validator


class Model(BaseModel):
    model_config = ConfigDict(extra='forbid', frozen=True)


class PathsConfig(Model):
    root_dir: DirectoryPath
    data_dir: Path
    output_dir: Path
    multirun_dir: Path


class SystemConfig(Model):
    inputfile: FilePath
    ions: Dict[str, int]
    reference_index: NonNegativeInt
    N_atom: PositiveInt
    # not FilePath: without them candidates are generated and nothing is excluded
    compositionfile: Optional[str] = None
    excludefile: Optional[str] = None
    limits: Optional[Dict[str, Tuple[NonNegativeInt, NonNegativeInt]]] = None

    @model_validator(mode='after')
    def limits_match_ions(self) -> 'SystemConfig':
        if self.limits is not None and set(self.limits) != set(self.ions):
            raise ValueError(f'limits must be given for exactly the ions {sorted(self.ions)}, got {sorted(self.limits)}')
        return self


class PathMode(Model):
    name: Literal['path']
    seeds_type: Literal['random', 'segmented']
    n_seeds: PositiveInt
    disect: PositiveInt
    max_iter: PositiveInt


class SuggestMode(Model):
    name: Literal['suggest']


class GenerateMode(Model):
    name: Literal['generate']


class BOConfig(Model):
    acquisition: Literal['qlogei', 'ts']
    kernel: Literal['matern', 'rbf']
    batch_size: PositiveInt


class PhaseBOConfig(Model):
    """Defaults live in config/; this only checks the composed configuration."""
    paths: PathsConfig
    system: SystemConfig
    mode: Annotated[Union[PathMode, SuggestMode, GenerateMode], Field(discriminator='name')]
    bo: BOConfig
    seed: Optional[int]
    show_plots: bool

    @classmethod
    def from_hydra(cls, cfg: DictConfig) -> 'PhaseBOConfig':
        return cls.model_validate(OmegaConf.to_container(cfg, resolve=True))

# Accelerated discovery of stable compositions in inorganic materials with Bayesian optimisation of the chemical phase fields

3.02.2021 Andrij Vasylenko

## Functionality

The modes of running the code are:

1) mode = 'path' can calculate a 'would-be' Bayesian Optimisation path
for the previously calculated compositions in the phase field

2) mode = 'suggest' suggests new compositions for calculations,
based on the precomputed results.

3) mode = 'generate' generates the list of candidate compositions that can be edited and further used
for limiting the candidates to the particular compositions only in a subsequent run in 'suggest' mode

## Requirements

Python 3.11–3.13 and [uv](https://docs.astral.sh/uv/).

## Dependencies:
numpy;
scipy;
pandas;
matplotlib;
pymatgen;
BoTorch (PyTorch);
Hydra

Dependencies are pinned in `uv.lock` and installed automatically.

## Installation
`uv sync`

## Usage
1) Prepare a 2-column table, where each row has a composition 
and its value of Total Energy as a .csv file.
Make sure you include reference compositions in the phase field.
2) Describe the phase field in `config/system/` (see `config/system/LiSnSCl.yaml`),
providing names of the files, atoms, and their oxidation states.
3) From anywhere in the repository, run:

`uv run phasebo`

and choose options or override any value of the configuration on the command line:

`uv run phasebo system=LiSnSCl mode=path mode.n_seeds=10 bo.acquisition=ts`

`uv run phasebo --help` lists the options and shows the configuration.

Each run writes its log (`phasebo.log`), results, plots and the configuration it used (`.hydra/`)
to `outputs/<system>/<date>/<time>_<mode>/` in the repository (`paths.output_dir` to change it).

To sweep over values, with one directory per run under `multirun/`:

`uv run phasebo -m mode=path bo.acquisition=qlogei,ts seed=0,1,2`

or, as a preset over 10 seeds, `uv run phasebo experiment=compare_acquisitions`.
Each 'path' run also writes `bo_path.csv`, and after a sweep `summary.csv` in the sweep directory
lists, for every 'path' run, the swept values, the number of stable compositions found
and the best and mean energy of the BO evaluations.

## Example
The default run, `uv run phasebo`, results in the outputs in `example`

## Reference
Please consider citing this tool:
A. Vasylenko et al,
Inferring energy–composition relationships with Bayesian optimization enhances exploration of inorganic materials.
The Journal of Chemical Physics 160, 5, 054110 (2024)


Bayesian optimisation is implemented with BoTorch
(earlier versions, including the paper above, used GPyOpt).
@inproceedings{balandat2020botorch,
  author =    {Balandat, Maximilian and Karrer, Brian and Jiang, Daniel R. and Daulton, Samuel and Letham, Benjamin and Wilson, Andrew Gordon and Bakshy, Eytan},
  title =     {{BoTorch: A Framework for Efficient Monte-Carlo Bayesian Optimization}},
  booktitle = {Advances in Neural Information Processing Systems 33},
  year =      {2020}
}

## Configuration

The configuration is composed from the files in `config/`: `config.yaml` and one option of each group,
`paths`, `system` (the phase field) and `mode`.
Paths start from the repository root, which is marked by the `.project-root` file.
Unknown keys and invalid values stop the run before it starts (see `src/phasebo/schema.py`).

 parameter | value 
---|--- 
*paths.data_dir*     | (default: data/ in the repository) Folder of the input files.
*paths.output_dir*   | (default: outputs/ in the repository) Where single runs are written.
*paths.multirun_dir* | (default: multirun/ in the repository) Where sweeps are written.
*system*       | (default: LiSnSCl) Phase field, a file in `config/system/`.
*system.inputfile*    | Input file. A table of compositions and their total energies.
*system.compositionfile*  | Input file. A list of candidate compositions (formulas) to consider. If not found, the candidates will be generated automatically.
*system.excludefile*  | Input file. A list of compostions (formulas) to exclude from convex hull calculations as well as from candidates. If not found, no candidates are excluded.
*system.reference_index* | Row of the inputfile where the reference compositions start.
*system.ions*         | Ions and oxidation states, e.g. {Li: 1, Sn: 4, S: -2, Cl: -1}.
*system.limits*       | Range of the number of atoms of each element in generated candidates.
*system.N_atom*       | Maximum number of atoms per unit cell in suggested compositions (in 'suggest' and 'generate' modes)
*mode*         | (default: suggest) Mode of calculations: the best path so far ('path'); suggest next compositions for CSP based on the available results ('suggest'); generate candidate compositions into candidates_list.csv ('generate') 
*mode.seeds_type*   | (default: 'random') Method to choose seeds in mode == 'path': 'segmented': Seeds are picked from a segmented phase field. 'random' seeds are selected randomly. 
*mode.n_seeds*      | (default: 23) Number of 'random' seeds in mode == 'path'.
*mode.disect*       | (default: 3) Number of sections of the phase field (disect x disect), from which the 'segmented' seeds are selected.
*mode.max_iter*     | (default: 10) Maximum number of iterations in mode == 'path'. 
*bo.batch_size*   | (default: 4) Number of compositions suggested per iteration.
*bo.acquisition*  | (default: 'qlogei') Batch selection: 'qlogei' (batch log expected improvement) or 'ts' (Thompson sampling, as in the GPyOpt version).
*bo.kernel*       | (default: 'matern') Gaussian process kernel: 'matern' (Matérn 5/2, as in the paper and the GPyOpt version) or 'rbf' (BoTorch's default).
*seed*         | (default: null) Random seed, for reproducible runs.
*show_plots*   | (default: true) Open plot windows in single runs; plots are always saved.

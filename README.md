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

python-3.9

## Dependencies:
numpy;
pandas;
pymatgen;
scikit-learn
GpyOpt

Dependencies can be installed automatically during installation.

## Installation
`pip install .`

## Usage
1) Prepare a 2-column table, where each row has a composition 
and its value of Total Energy as a .csv file.
Make sure you include reference compositions in the phase field.
2) Modify input_config.yaml file accordingly, 
providing names of the file, atoms, and their oxidation states.
3) run:

`python -m phasebo`

will use input_config.yaml by default

or 

`python -m phasebo --config path/to/my_config.yaml`

## Example
The default run with input_config.yaml results int the outputs in `example`

## Reference
Please consider citing this tool:
"A. Vasylenko et al 


The code is based on GPyOpt implementation of Bayesian optimisation 
algorithm.
@Misc{gpyopt2016,
  author =   {The GPyOpt authors},
  title =    {GPyOpt: A Bayesian Optimization framework in Python},
  howpublished = {\url{http://github.com/SheffieldML/GPyOpt}},
  year = {2016}
}

## Parameters of the input configuration file 

 parameter | value 
---|--- 
*mode*         | (default: 'suggest') Mode of calculations: the best path so far ('path'); suggest next compositions for CSP based on the available results ('suggest'); generate candidate compositions into candidates_list.csv ('generate') 
*inputfile*    | (default: LiSnSCl_700eV.csv) Input file. A table of compositions and their total energies.
*compositionfile*  | (default: None) Input file. A list of candidate compositions (formulas) to consider. If not provided, the candidates will be generated automatically.
*excludefile*  | (default: None) Input file. A list of compostions (formulas) to exclude from convex hull calculations as well as from candidates. If not provided, no candidates are excluded.
*ions*         | (default: {'Li':1,'Sn':4,'S':-2,'Cl':-1}) Ions and oxidation states.
*seeds_type*   | (default: 'random') Method to choose seeds in mode == 'path': 'segmented': Seeds are picked from a segmented phase field. 'random' seeds are selected randomly. 
*disect*       | (default: 4) Number of sections of the phase field (disect x disect), from which the 'segmented' seeds are selected.
*N_atom*       | (default: 24) Maximum number of atoms per unit cell in suggested compositions (in 'suggest' and 'generate' modes)
*max_iter*     | (default: 10) Maximum number of iterations. 

# Exploration of the phase fields' minima with Bayesian Optimisation

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

python-3.8

pip (version 19.0 or later)

OS:

Ubuntu (version 18.04 or later)

MacOS (Catalina 10.15.6 or later)

## Dependencies:
numpy;
pandas;
pymatget;
GpyOpt

Dependencies can be installed automatically during installation.

## Installation
`pip install phasebo`
or from the source:
`pip install .`

## Usage
1) Prepare a 2-column table, where each row has a composition 
and its value of Total Energy as a .csv file.
Make sure you include reference compositions in the phase field.
2) Modify input.py file accordingly, 
providing names of the file, atoms, and their oxidation states.
3) run:

`python input.py`



## Reference
The code is based on GPyOpt implementation of Bayesian optimisation 
algorithm.
@Misc{gpyopt2016,
  author =   {The GPyOpt authors},
  title =    {GPyOpt: A Bayesian Optimization framework in Python},
  howpublished = {\url{http://github.com/SheffieldML/GPyOpt}},
  year = {2016}
}

## Parameters of the input file 
(default values are in the __main__.py file)

 parameter | value 
---|--- 
 *mode*         | (default: 'path') Mode of calculations: the best path so far ('path'); suggest next compositions for CSP based on the available results ('suggest'); generate candidate compositions into candidates_list.csv ('generate') 
*ifile*  | (default: LiSnSCl_700eV.csv) Input file. A table of compostions and their total energies.
*cfile*  | (default: None) Input file. A list of candidate compostions (formulas) to consider. If not provided, the candidates will be generated automatically.
*ions*   | (default: {'Li':1,'Sn':4,'S':-2,'Cl':-1}) Ions and oxidation states.
*seeds*  | (default: 'segmented') Method to choose seeds in mode == 'path': 'segmented': Seeds are picked from a segmented phase field. 'random' seeds are selected randomly - decreased efficiency. 
*disect* | (default: 4) Number of sections of the phase field (disect x disect), from which the 'segmented' seeds are selected.
*N_atom* | (default: 24) Maximum number of atoms per unit cell in suggested compositions (in 'suggest' and 'generate' modes)
*max_iter* | (default: 10) Maximum number of iterations. 

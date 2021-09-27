from phasebo.__main__ import run
import pandas as pd

#------------- INPUT ---------------

ifile = 'test_data.csv'               # File with compositions and total enthalpies
                                      # 
ions = {'Li':1,'Mg':2,'Al':3,'P':5,'O':-2} # Ions and oxidation states
                                      # 
mode = 'path'                         # Modes for running BO:
                                      #    'path':     Calculates a would-be-BO-path towards 
                                      #                a composition with minimum E above convex hull
                                      #    'suggest':  Calculates next best suggested compositions 
                                      #    'generate': Generates a list of candidate compositions
                                      # 
seeds = 'random'                      # Method to choose seeds in mode == 'path':
                                      #    'random' seeds are selected randomly
                                      #    'segmented': (deprecated) seeds are picked from sections of the phase field.
                                      # 
n_seeds = 23                           # Number of seeds. default: 9 
                                      # 
limits = None
#{'Li': [1, 10],              # Min and max limits on amounts of atoms of each atomic element 
#          'Sn': [1, 5],               # the limits will be used when candidates are generated in 'suggest' and 'generate' 
#           'S': [1, 10],              # modes.
#          'Cl': [1, 5]                # If limits are not specified (default), all charge-balanced compositions will be generated
#         }                            #

N_atom = 24                           # Maximum number of atoms per unit cell in suggested compositions
max_iter = 10                         # Evaluation budget for BO
                                      # 
                                      # 
#------------- READ ----------------- #
log = "logfile"                       # Name of the file, where results will be printed with a timestamp
df = pd.read_csv(ifile, header=0)     #  
compositions = df.values              # Select candidate and reported compositions
references = df.values[:112]          #
                                      # Read candidate compositions file if provided
try:
    next_formulas = [i[0] for i in pd.read_csv(cfile).values]
except Exception:
    next_formulas = None

#------------- PLOT ----------------
plot_mode = 'screen'                  # Select plotting mode: 'web' creates interactive plot in a browser (deprecated)
                                      # 'screen' (default) plots on a screen. 

#------------- RUN -----------------  # You shouldn't need to modify run function call or anything in phase_fields_bo/ folder

run(compositions, references, ions, mode, N_atom, seeds, n_seeds, max_iter, log, limits, next_formulas=next_formulas, allow_negative=False)

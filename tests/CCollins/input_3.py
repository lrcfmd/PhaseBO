from phasebo.__main__ import run
import pandas as pd

#------------- INPUT ---------------

ifile = 'bo_ref_data.csv'           # File with compositions and total enthalpies
#cfile = 'candidates_list.csv'         # File with candidate compositions;
                                      # if not provided (default), the candidates will be generated automatically
                                      # 
ions = {'Li':1,'B':3,'Zn':2,'O':-2,'S':-2} # Ions and oxidation states
                                      # 
mode = 'generate'                         # Modes for running BO:
                                      #    'path':     Calculates a would-be-BO-path towards 
                                      #                a composition with minimum E above convex hull
                                      #    'suggest':  Calculates next best suggested compositions 
                                      #    'generate': Generates a list of candidate compositions
                                      # 
seeds = 'random'                      # Method to choose seeds in mode == 'path':
                                      #    'random' seeds are selected randomly
                                      #    'segmented': (deprecated) seeds are picked from sections of the phase field.
                                      # 
n_seeds = 100                           # Number of seeds. default: 9 
                                      # 
limits = {'Li': [0, 5],              # Min and max limits on amounts of atoms of each atomic element 
          'B':  [0, 5],               # the limits will be used when candidates are generated in 'suggest' and 'generate' 
          'Zn': [0, 5],              # modes.
          'O':  [0, 5],
	  'S':  [0, 5]				# If limits are not specified (default), all charge-balanced compositions will be generated
         }                            #

#limits=None

N_atom = 35                           # Maximum number of atoms per unit cell in suggested compositions
max_iter = 100                         # Evaluation budget for BO
                                      # 
                                      # 
#------------- READ ----------------- #
log = "logfile"                       # Name of the file, where results will be printed with a timestamp
df = pd.read_csv(ifile, header=0)     #  
compositions = df.values              # Select candidate and reported compositions
references = df.values[195:]          #
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


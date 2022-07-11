from phase_fields_bo.__main__ import run
import pandas as pd

#------------- INPUT ---------------

ifile = 'LiSiSF_xtalopt.csv'         # File with compositions and total enthalpies
ions = {'Li':1,'Si':4,'S':-2,'F':-1} # Ions and oxidation states
                                      # 
#mode = 'suggest'                         # Modes for running BO:
mode = 'generate'                                      #    'path':    Calculates a would-be-BO-path towards 
                                      #               a composition with minimum E above convex hull
                                      #    'suggest': Calculates next best suggested compositions 
                                      # 
seeds = 'random'                      # Method to choose seeds in mode == 'path':
                                      #    'random' seeds are selected randomly
                                      #    'segmented': (deprecated) seeds are picked from sections of the phase field.
                                      # 
N_atom = 40                           # Maximum number of atoms per unit cell in suggested compositions
max_iter = 10                         # Evaluation budget for BO

#------------- READ -----------------
log = "logfile"                       # Name of the file, where results will be printed with a timestamp
df = pd.read_csv(ifile, header=0)     #  
compositions = df.values              # Select candidate and reported compositions
references = df.values[22:]          #

#------------- PLOT ----------------
plot_mode = 'screen'                  # Select plotting mode: 'web' creates interactive plot in a browser
                                      # 'screen' (default) plots on a screen. 

#------------- RUN -----------------  # You shouldn't need to modify run or anything in phase_fields_bo/ folder

run(compositions, references, ions, mode, N_atom, seeds, max_iter, log, allow_negative=False)

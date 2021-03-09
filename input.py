from phase_fields_bo.__main__ import run
import pandas as pd

#------------- INPUT ---------------

ifile = 'LiSnSCl_700eV.csv'           # File with compositions and total enthalpies
ions = {'Li':1,'Sn':4,'S':-2,'Cl':-1} # Ions and oxidation states
                                      # 
mode = 'path'                         # Modes for running BO:
                                      #    'path':    Calculates a would-be-BO-path towards 
                                      #               a composition with minimum E above convex hull
                                      #    'suggest': Calculates next best suggested compositions 
                                      # 
seeds = 'random'                      # Method to choose seeds in mode == 'path':
                                      #    'random' seeds are selected randomly
                                      #    'segmented': (redundant) Seeds are picked from sections of the phase field.
                                      # 
N_atom = 24                           # Maximum number of atoms per unit cell in suggested compositions
max_iter = 10                         # Evaluation budget for BO

#------------- READ -----------------
log = "logfile"                       # Name of the file, where results will be printed 
df = pd.read_csv(ifile, header=0)     #  
compositions = df.values              # Select candidate and reported compositions
references = df.values[195:]          #

#------------- PLOT ----------------
plot_mode = 'screen'                  # Select plotting mode: 'web' creates interactive plot in a browser
                                      # 'screen' (default) plots on a screen. 

#------------- RUN -----------------  # You shouldn't need to modify run or anything in phase_fields_bo/ folder

run(compositions, references, ions, mode, N_atom, seeds, disect, max_iter, log)

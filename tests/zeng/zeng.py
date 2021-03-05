from phase_fields_bo.__main__ import run
import pandas as pd

#------------- INPUT ---------------

ifile = 'bayesian.csv'           # File with compositions and total enthalpies
ions = {'Ba':2,'Nb':4,'Mg':2,'O':-2}
                                      # 
mode = 'suggest'                         # Modes for running BO:
                                      #    'path':    Calculates a would-be-BO-path towards 
                                      #               a composition with minimum E above convex hull
                                      #    'suggest': Calculates next best suggested compositions 
                                      # 
seeds = 'random' #segmented'                   # Method to choose seeds in mode == 'path':
                                      #    'segmented': Seeds are picked from a segmented phase field.
                                      #    'random' seeds are selected randomly - decreased efficiency 
                                      # 
disect = 3                            # Number of sections of the phase field: disect x disect 
N_atom = 24                           # Maximum number of atoms per unit cell in suggested compositions
max_iter = 2                         # Evaluation budget for BO

#------------- READ -----------------
log = open("logfile", 'a')            # 
df = pd.read_csv(ifile, header=0)     #  
compositions = df.values              # Select candidate and reported compositions
references = df.values[21:]          #


#------------- RUN -----------------
run(compositions, references, ions, mode, N_atom, seeds, disect, max_iter, log)

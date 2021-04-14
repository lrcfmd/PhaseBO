import sys
import time
import numpy as np
import pandas as pd
from GPyOpt.methods import BayesianOptimization
from phase_fields_bo.phase_field_bo import PhaseFieldBO

def run(compositions, references, ions, mode, Ntot, seeds, max_iter, log, next_formulas=None, allow_negative=False):
    """ BO run """
    bopt = PhaseFieldBO(compositions,
                    references, 
                    ions, 
                    mode, 
                    seeds,
                    exclude_zeros=True,
                    Ntot=Ntot,
                    max_iter=max_iter,
                    next_formulas=next_formulas,
                    allow_negative=allow_negative)

    # --------- PLOT & PRINT -------------

    convex = bopt.plot_convex()
    convex.show()

    # --------- logtime
    t = time.localtime()
    timestamp = time.strftime('%b-%d-%Y_%H%M', t)
    log = open(f"{log}-{timestamp}",'a')

    if mode == 'path':
       # bopt.plot_path() -- to add {pd_coords: 2d_square_coords} for plotting on 2d square
        bopt.bo.plot_convergence()
        bopt.bo.plot_acquisition()
        bopt.print_results(log)
    elif mode == 'suggest':
       # bopt.plot_suggested() -- to add {pd_coords: 2d_square_coords} for plotting on 2d square
        bopt.print_results(log)

    return bopt

if __name__ == '__main__':
    
    ifile = 'LiSnSCl_700eV.csv'           # File with compositions and total enthalpies
    ions = {'Li':1,'Sn':4,'S':-2,'Cl':-1} # Ions and oxidation states
                                          # 
    mode = 'suggest'                         # Modes for running BO:
                                          #    'path':    Calculates a would-be-BO-path towards 
                                          #               a composition with minimum E above convex hull
                                          #    'suggest': Calculates next best suggested compositions 
                                          # 
    seeds = 'random'                      # Method to choose seeds in mode == 'path':
                                          #    'random' seeds are selected randomly 
                                          #    'segmented': Seeds are picked from a segmented phase field.
                                          #            Supports only phase fields represented as a square.
                                          # 
    N_atoms = 24                          # Maximum number of atoms in a unit cell in suggested compositions
    max_iter = 10                         # Evaluation budget for BO
    log = "logfile"                       # 
    df = pd.read_csv(ifile, header=0)     #  
    compositions = df.values              # Select candidate and reported compositions
    references = df.values[195:]          #
    candidatesf = 'candidates_list.csv'
    next_formulas = [i[0].strip() for i in pd.read_csv(candidatesf).values]

    run(compositions, references, ions, mode, N_atoms, seeds, max_iter, log, next_formulas=None, allow_negative=False)

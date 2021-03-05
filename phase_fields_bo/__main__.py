import time
import numpy as np
import pandas as pd
from GPyOpt.methods import BayesianOptimization
from phase_fields_bo.phase_field_bo import PhaseFieldBO

def run(compositions, references, ions, mode, Ntot, seeds, disect, max_iter, log):
    """ BO run """
    bopt = PhaseFieldBO(compositions,
                    references, 
                    ions, 
                    mode, 
                    seeds,
                    exclude_zeros=True,
                    n_seeds=9,
                    disect=disect,
                    Ntot=Ntot,
                    max_iter=max_iter)

    # --------- PLOT & PRINT -------------

    convex = bopt.plot_convex()
    convex.show()

    if mode == 'path':
        bopt.plot_path()
        bopt.bo.plot_convergence()
        bopt.bo.plot_acquisition()
        bopt.print_results(log)
    elif mode == 'suggest':
        bopt.plot_suggested()
        t = time.localtime()
        timestamp = time.strftime('%b-%d-%Y_%H%M', t)
        log = open(f"{log}-{timestamp}",'a')
        bopt.print_results(log)

if __name__ == '__main__':
    
    ifile = 'LiSnSCl_700eV.csv'           # File with compositions and total enthalpies
    ions = {'Li':1,'Sn':4,'S':-2,'Cl':-1} # Ions and oxidation states
                                          # 
    mode = 'suggest'                         # Modes for running BO:
                                          #    'path':    Calculates a would-be-BO-path towards 
                                          #               a composition with minimum E above convex hull
                                          #    'suggest': Calculates next best suggested compositions 
                                          # 
    seeds = 'segmented'                   # Method to choose seeds in mode == 'path':
                                          #    'segmented': Seeds are picked from a segmented phase field.
                                          #    'random' seeds are selected randomly - decreased efficiency 
                                          # 
    disect = 4                            # Number of sections of the phase field: disect x disect 
    N_atoms = 24                          # Maximum number of atoms in a unit cell in suggested compositions
    max_iter = 10                         # Evaluation budget for BO
    log = "logfile"                       # 
    df = pd.read_csv(ifile, header=0)     #  
    compositions = df.values              # Select candidate and reported compositions
    references = df.values[195:]          #


    run(compositions, references, ions, mode, N_atoms, seeds, disect, max_iter, log)

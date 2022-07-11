from phasebo.__main__ import run
import pandas as pd

#------------- INPUT ---------------

ifile = 'LiBSI.csv'
#cfile = 'candidates_list.csv'
ions = {'Li':1,'B':3,'I':-1,'S':-2} # Ions and oxidation states
                                      # 
mode = 'generate'                         # Modes for running BO:
#mode = 'suggest'                         # Modes for running BO:
                                      #    'path':    Calculates a would-be-BO-path towards 
                                      #               a composition with minimum E above convex hull
                                      #    'suggest': Calculates next best suggested compositions 
                                      # 
#limits = {'Li':[1,20],
#          'B': [1,4],
#          'S': [1,10],
#          'I': [1,20]}

seeds = 'random'                      # Method to choose seeds in mode == 'path':
                                      #    'random' seeds are selected randomly
                                      #    'segmented': (deprecated) seeds are picked from sections of the phase field.
                                      # 
N_atom = 40                           # Maximum number of atoms per unit cell in suggested compositions
max_iter = 10                         # Evaluation budget for BO

#------------- READ -----------------
log = "logfile"                       # Name of the file, where results will be printed with a timestamp
df = pd.read_csv(ifile, header=0, comment="#")     #  
compositions = df.values              # Select candidate and reported compositions
references = df.values[:19]          #
try:
    next_formulas = [i[0] for i in pd.read_csv(cfile).values]
except Exception:
    next_formulas = None


#------------- PLOT ----------------
plot_mode = 'screen'                  # Select plotting mode: 'web' creates interactive plot in a browser
                                      # 'screen' (default) plots on a screen. 


#------------- RUN -----------------  # You shouldn't need to modify run or anything in phase_fields_bo/ folder

bopt = run(compositions, references, ions, mode, N_atom, seeds, 9,  max_iter, log, limits=limits, next_formulas=None, allow_negative=False)

#------------- POSTPROCES -----------
#for entry in bopt.computed_entries:
#    print(entry.composition, entry.energy, 1000*bopt.pd.get_decomp_and_e_above_hull(entry,allow_negative=True)[1])
#    print(entry.composition, entry.energy, 1000*bopt.pd.get_e_above_hull(entry))

#  if 1000*bopt.pd.get_e_above_hull(entry) == 0:
#      print(entry.composition, entry.energy, 1000*bopt.pd.get_equilibrium_reaction_energy(entry))
#  else:
#      print(entry.composition, entry.energy, 1000*bopt.pd.get_e_above_hull(entry))


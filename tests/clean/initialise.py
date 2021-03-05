import numpy as np
from GPyOpt.methods import BayesianOptimization
from convexhull import get_convexhull
import list_compositions
import matplotlib 
from plotter import get_plot
import sys
#------ INPUT --------------------------------

Ntot = 20                             # number of atoms in UC to consider
ions = {'Li':1,'Zn':2,'S':-2,'Cl':-1} # ions and oxidation states
references_list = 'reference.dat'     # list of reference compositions with energies
compositions_list = "list3.dat"       # list of compositions with energies
#log = open("results_and_next", 'a')   # file that will contain the results
log = "results_and_next_1"   # file that will contain the results
batch = 35                            # a batch size to choose from in Thompson sampling

# --------------------------------------------


elements = list(ions.keys())

# get initial compositions, energies and elements from the compositions_list
compositions, energies = list_compositions.initial(compositions_list)

# plot initial PES, and get stoichiometries(var1, var2) and PES energies, and phase diagram
#var1, var2, pes = get_convexhull(references_list, compositions, energies, elements, log, allref=True)
var1, var2, pes = get_convexhull(references_list, compositions, energies, elements, allref=True)
test = np.array([var1, var2]).T
print(test.shape, pes.shape)
sys.exit(0)
# create domain
dom, candidates = list_compositions.generate(ions, compositions, Ntot)
domain = [{'name': 'var_1', 'type': 'bandit', 'domain': dom}]

# run optimization
bo_step = BayesianOptimization(f = None, 
                               domain = domain,
                               X = np.array([var1, var2]).T,
                               Y = pes.reshape(-1,1),
                               evaluator_type = 'thompson_sampling',
                               batch_size = batch)

x_next = bo_step.suggest_next_locations()

log = list_compositions.print_pes(compositions, energies, log)
list_compositions.print_next(x_next, candidates, log)

get_plot(elements, np.array([var1, var2, pes]).T, x_next, dom)

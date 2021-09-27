import sys
import numpy as np
import random
from itertools import product
from numpy.random import seed
from GPyOpt.methods import BayesianOptimization
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity as similarity
from sklearn.preprocessing import StandardScaler as SS

NSEEDS = 16
NBATCH = 10
NRUNS = 50

df = pd.read_csv('Li-Mg-Al-P-O.csv', comment='#')

# variables
Li=df['f(Li)'].values
Mg=df['f(Mg)'].values
Al=df['f(Al)'].values
P= df['f(P)'].values
O= df['f(O)'].values

references = np.array([np.array([a,b,c,d,e]) for a, b, c, d, e in zip(Li, Mg, Al, P, O)])

# result values
obj = df['E (meV/atom)'].values
obj = obj.reshape(-1,1)
# =================================   

def f(x):
    """ function of ehull energies for all fractional coordinates """
    return abs(np.dot(np.where(np.all(references == x, axis=1), 1, 0), obj))

def l2s(l):
    return ','.join(map(str, [round(i,3) for i in l]))
# =================================   

# 10 independent runs
for r in range(1,11):

    # random seeds
    rseeds = np.random.choice(np.arange(len(references)), NSEEDS)
    seeds = references[rseeds]
    seeds_s = [l2s(s) for s in seeds]

    # candidates - observed experiments
    candidates = np.array([r for r in references if l2s(r) not in seeds_s])

    X_init = np.array([np.array(s) for s in seeds])
    Y_init = obj[rseeds]

    domain = [{'name': 'var_1', 'type': 'bandit', 'domain': candidates}]

    # Run optimisation
    bo = BayesianOptimization(f=f,
               domain=domain,
               X = X_init,
               Y = Y_init,
               evaluator_type = 'thompson_sampling',
               batch_size = NBATCH,
               de_duplication = True)

    bo.run_optimization(NRUNS, verbosity=False)  # find a path in a known domain
    #bo.plot_convergence()

    print ('RUN:', r, file=open('logfile', 'a'))
    print ('SEEDS:', file=open('logfile', 'a'))
    for s, y in zip(X_init, Y_init):
        print(s, y[0], file=open('logfile', 'a'))

    results = np.array([f(x)[0] for x in bo.X])
    print ('BEST:', min(results), 'ITERATIONS:', 1+np.argmin(results), file=open('logfile', 'a')) 
    print(50*'=', file=open('logfile', 'a'))

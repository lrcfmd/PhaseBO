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

# 10 independent runs
AR = 0
for r in range(1,11):

    e = 0; i=0
    while e != 17.5: 
        e = np.random.choice(obj)
        i += 1
    print ('RANDOM choice:', i)
    AR += i

print ('Average:', AR/10)


import os, sys
import pandas as pd
import numpy as np	
from itertools import product as P
from pprint import pprint
#from iteration_utilities import deepflatten
from pymatgen.entries.computed_entries import ComputedEntry

def initial(compositions_list):
    c_list = open(compositions_list).readlines()
    compositions = [i.split()[0] for i in c_list]
    energies = [float(i.split()[1]) for i in c_list]
    return compositions, energies

def span(ions, Ntot, limits=None):
    if limits:
        amounts = [np.arange(limits[i][0], limits[i][1]) for i in ions]
    else:
        amounts = [np.arange(1,Ntot+1) for i in ions]
    PP = list(P(*amounts))    
    return [i for i in PP if sum(i) <= Ntot]  

def balance(amounts, ions):
    balanced = []
    for a in amounts:
        if sum(a) == 0: continue 
        if sum([n * charge for n,charge in zip(a, ions)]) == 0:
            balanced.append(a)
    return balanced

def generate(ions, inlist, Ntot, limits):
    amounts = span(ions, Ntot, limits)
    amounts = balance(amounts, list(ions.values()))
    symbols = list(ions.keys())
    names = []  

    for amount in amounts: 
        name = []
        for s, n in zip(symbols, amount):
            name.append(f'{s}{n}')

        if ''.join(name) not in inlist:
            names.append(''.join(name))

   
    return names

def print_pes(compositions, energies, log):
    if os.path.exists(log):
        if log[-1].isdigit(): 
            i  = int(log[-1]) + 1
            log = log[:-1] + str(i)
        else: log += '0'

    for c,e in zip(compositions, energies):
        print(c, e, file=open(log, 'a'))

    return log 

def print_next(x_next, candidates, log):
    x_next = [' '.join(map(str,x)) for x in list(x_next)]
    for x in set(x_next):
        next_formula = candidates[x]
        print('next:', next_formula, file=open(log, 'a'))

if __name__=="__main__":
    ions = {'Li':1, 'Zn':2, 'S':-2, 'Cl':-1}
    #ions = {'Ba':2,'Nb':4,'Mg':2,'O':-2}
    Ntot = 24

    compositions = []
    names = generate_recursive(ions, compositions, Ntot)
    pprint(names)

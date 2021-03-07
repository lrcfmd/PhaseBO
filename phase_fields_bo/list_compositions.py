import os
import pandas as pd
import numpy as np	
from itertools import product as P
from pymatgen.entries.computed_entries import ComputedEntry

def initial(compositions_list):
    c_list = open(compositions_list).readlines()
    compositions = [i.split()[0] for i in c_list]
    energies = [float(i.split()[1]) for i in c_list]
    return compositions, energies

def Span(ions, Ntot):
    amounts = [np.arange(1,Ntot+1) for i in ions]
    def PP(n):
        if n == 0:
            return P(amounts[n])
        else:
            return P(amounts[n], PP(n-1))

    def unnest(amounts):
        result = []
        for i in amounts:
            if isinstance(i,tuple):
                 result.extend(unnest(i))
            else:
                 result.append(i)
        return result

    return [unnest(i) for i in PP(len(ions)-1) if sum(unnest(i)) <= Ntot]

def balance(amounts, ions):
    balanced = []
    for a in amounts:
        if sum([n * charge for n,charge in zip(a, ions)]) == 0:
            balanced.append(a)
    return balanced

def generate(ions, inlist, Ntot):
    ''' Lists all charge-ballanced quaternary compositions with Natoms < Ntot and not in inlist '''
    print (f'Generating candidate compositions with Natoms <= {Ntot}')
    el = [i for i in ions]
    newlist = {}
    # discrete variables (number of atoms of species) for BO
    var1, var2 = [], []

    for a1 in range(1, Ntot):
        for a2 in range(1, Ntot):
            for c1 in range(1, Ntot):
                for c2 in range(1, Ntot):
                    charges = [n * a for n,a in zip([a1,a2,c1,c2], list(ions.values()))]
                    if sum(charges) or a1 + a2 + c1 + c2 > Ntot: continue

                    name = f"{el[0]}{a1}{el[1]}{a2}{el[2]}{c1}{el[3]}{c2}"
                    
                    if ComputedEntry(name, 0).composition not in inlist:
                        x = round(c1/(c1+c2), 3)
                        y = round(a1/(a1+a2), 3)
                        var1.append(x)
                        var2.append(y)
                        newlist[f'{x} {y}'] = name
    return np.array([tuple([i,j]) for i,j in zip(var1,var2)]), newlist


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

def unnest(amounts):
    result = []
    for i in amounts:
        if isinstance(i,tuple):
             result.extend(unnest(i))
        else:
             result.append(i)
    return result


if __name__=="__main__":
    import time
    ions = {'Li':1, 'Zn':2, 'S':-2,'Cl':-1}
    #ions = {'Ba':2,'Nb':4,'Mg':2,'O':-2}
    Ntot = 24
    compositions = []
    s = time.time()
    test1, test2 = generate(ions, compositions, Ntot)
    print('Loop time:', time.time()-s)

    s = time.time()
    amounts = Span(ions, Ntot)
    amounts = balance(amounts, list(ions.values()))
    print('Recursion time:', time.time()-s)

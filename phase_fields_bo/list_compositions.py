import os
from numpy import array
from pymatgen.entries.computed_entries import ComputedEntry

def initial(compositions_list):
    c_list = open(compositions_list).readlines()
    compositions = [i.split()[0] for i in c_list]
    energies = [float(i.split()[1]) for i in c_list]
    return compositions, energies


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
                    
                    print(name, charges)
                    if ComputedEntry(name, 0).composition not in inlist:
                        x = round(c1/(c1+c2), 3)
                        y = round(a1/(a1+a2), 3)
                        var1.append(x)
                        var2.append(y)
                        newlist[f'{x} {y}'] = name
    return array([tuple([i,j]) for i,j in zip(var1,var2)]), newlist


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
    ions = {'Li':1,'Zn':2,'S':-2,'Cl':-1}
    #ions = {'Ba':2,'Nb':4,'Mg':2,'O':-2}
    Ntot = 10
    compositions = []
    test1, test2 = generate(ions, compositions, Ntot)

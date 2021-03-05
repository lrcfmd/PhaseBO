from numpy import array
import os

def initial(compositions_list):
    c_list = open(compositions_list).readlines()
    compositions = [i.split()[0] for i in c_list]
    energies = [float(i.split()[1]) for i in c_list]
    return compositions, energies


def generate(ions, inlist, Ntot):
    ''' Lists all charge-ballanced quaternary compositions with Natoms < Ntot and not in inlist '''
    print (f'Generating candidate compositions with Natoms <= {Ntot}')
    el = [i for i in ions]
    anions = [i for i in el if ions[i] < 0]
    cations = [i for i in el if ions[i] > 0] 
    newlist = {}
    # discrete variables (number of atoms of species) for BO
    var1, var2 = [], []
    for a1 in range(1, Ntot):
        for a2 in range(1, Ntot):
            nanions = sum([n * a for n,a in zip([a1, a2], [ions[i] for i in anions])])
            for c1 in range(1, Ntot):
                for c2 in range(1, Ntot):
                    ncations = sum([n * a for n,a in zip([c1, c2], [ions[i] for i in cations])])
                    if ncations != -nanions or a1 + a2 + c1 + c2 > Ntot: continue
                    name = f"{el[0]}{c1}{el[1]}{c2}{el[2]}{a1}{el[3]}{a2}"
                    if name not in inlist:
                        var1.append(c1/(c1+c2))
                        var2.append(a1/(a1+a2))
                        newlist[f'{var1[-1]} {var2[-1]}'] = name
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
    Ntot = 17
    inlist = open("list.dat").readlines()
    compositions = [i.split()[0] for i in inlist]

    _, __ = generate_atoms(ions, compositions, Ntot)

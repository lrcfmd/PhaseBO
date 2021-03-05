import numpy as np
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry

def get_reference(references_file):
    reffile = open(references_file, 'r').readlines()
    reference = []

    for line in reffile:
        c = line.split()[0]
        e = float(line.split()[1])
        reference.append(ComputedEntry(f'{c}', e))
    return reference

#def get_convexhull(references_file, compositions, energies, elements, log, allref=False):
def get_convexhull(references_file, compositions, energies, elements, allref=False):
    print('Calculating convex hull:')
    
    reference = get_reference(references_file)
    structures = []
    for c,e in zip(compositions, energies):
       structures.append(ComputedEntry(f'{c}', e))

    allstructures =  reference + structures
    PDfull=PhaseDiagram(allstructures)

    var1, var2, pes = [], [], []
    
    if allref: struc = allstructures
    else: struc = structures

    for entry in struc:
        c1 = entry.composition[elements[0]]
        c2 = entry.composition[elements[1]]
        a1 = entry.composition[elements[2]]
        a2 = entry.composition[elements[3]]
        if c1 + c2 != 0:
            var1.append(c1/(c1+c2))
        else: var1.append(1)
        if a1 + a2 != 0:
            var2.append(a1/(a1+a2))
        else: var2.append(1)
        pes.append(1000*PDfull.get_e_above_hull(entry))
#        print(str(entry.composition).replace(" ",""), 1000*PDfull.get_e_above_hull(entry), file=log)

    return np.array(var1), np.array(var2), np.array(pes)


# pullFromMP.py
 
import os.path
import sys
from pymatgen.ext.matproj import MPRester
import ase
 
mp = MPRester("4r0hw9Mm0wjGgkSeIE") # your mp key
phases = ['Li', 'B', 'S', 'Br']
 
##############################################
structs = mp.get_entries_in_chemsys(phases)
ids = [s.entry_id for s in structs]
 
for i in ids:
     n = 0
     entry = mp.get_entry_by_material_id(i, inc_structure='final')
     if 'Br' not in entry.composition: continue 

#    if 'Zr' not in entry.composition and 'Br' not in entry.composition:
#       print(entry.composition)
#       continue
     name = str(entry.composition).replace(" ", "") + '_POSCAR' + str(n)
     print(i, name)
     while os.path.isfile(name):
         n += 1
         name = name[:-1] + str(n)
     entry.structure.to(fmt="vasp", filename=name)
# end of script

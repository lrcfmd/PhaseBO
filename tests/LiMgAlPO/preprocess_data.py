import pandas as pd
import numpy as np

df = pd.read_csv('Li-Mg-Al-P-O.csv')

Li=df['Li'].values
Mg=df['Mg'].values
Al=df['Al'].values
P=df['P'].values
O=df['O'].values

elements=['Li','Mg','Al','P','O']
cs = np.asarray([Li,Mg,Al,P,O]).transpose()
formulas = []

for ck in cs:
    formula = ''
    for c,e in zip(ck, elements):
        if c != 0:
            formula += f'{e}{c}'
    print(formula)
    formulas.append(formula)

newd = pd.DataFrame({'Composition': formulas, 'Energy': df['E (meV/atom)'].values})
newd.to_csv('LiMgAlPO_forBO.csv', index=None)


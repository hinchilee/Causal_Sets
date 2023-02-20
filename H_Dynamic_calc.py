# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 20:07:26 2023

@author: leehi
"""

import pandas as pd
import itertools
from causalsetfunctions import n_sphere_surfacearea, find_entropy

T = 3
d_array = [2,3,4]
for d in d_array: 
    df = pd.read_csv(f'H_Dynamic{d}d.csv', names=['rho', 'H'], header=None)
    #print(df)
    rho_array = df['rho'].unique()
    rho_array.sort()
    for rho in rho_array:
        
        print('\n sprinkling density',rho, f'in {d} dimensions')
        iterationsNo = df[df['rho'] == rho].shape[0]
        print(f'number of iterations: {iterationsNo}')
        print(f'T = 3')
        totalHarray = [sum(x) for x in itertools.zip_longest(*[[int(x) if x != '' else 0 for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)]
        print('total Harray \n', totalHarray)
        
        totalLinks = 0 
        for i, value in enumerate(totalHarray):
            totalLinks += (i+1)*value
            
        print(f'Total Links: {totalLinks}')
        
        empiricalavalue = rho**((2-d)/d)*(totalLinks/iterationsNo)/ n_sphere_surfacearea(n = d - 2, r = T)
        
        print(f'Empirical a value {empiricalavalue}')
        
        entropy = find_entropy(totalHarray)
        print(f'Entropy: {entropy}')
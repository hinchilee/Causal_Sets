# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 20:07:26 2023

@author: leehi
"""
import numpy as np 
import pandas as pd
import itertools
from causalsetfunctions import n_sphere_surfacearea, find_entropy, linear
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

T = 3
d_array = [2,3,4]
for d in d_array: 
    df = pd.read_csv(f'H_Dynamic{d}d.csv', names=['rho', 'H'], header=None)
    #print(df)
    rho_array = df['rho'].unique()
    rho_array.sort()
    y_entropyList = list() 
    x = list()
    for rho in rho_array:
        
        print('\n sprinkling density',rho, f'in {d} dimensions')
        iterationsNo = df[df['rho'] == rho].shape[0]
        print(f'number of iterations: {iterationsNo}')
        print(f'T = 3')
        totalHarray = [sum(x) for x in itertools.zip_longest(*[[int(x) if x != '' else 0 for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)]
        print('total Harray \n', totalHarray)
        print('p_i:', totalHarray/ np.sum(totalHarray))
        
        totalLinks = 0 
        for i, value in enumerate(totalHarray):
            totalLinks += (i+1)*value
            
        print(f'Total Links: {totalLinks}')
        
        empiricalavalue = rho**((2-d)/d)*(totalLinks/iterationsNo)/ n_sphere_surfacearea(n = d - 2, r = T)
        
        print(f'Empirical a value {empiricalavalue}')
        
        entropy = find_entropy(totalHarray, iterationsNo)
        print(f'Entropy: {entropy}')
        
        y_entropyList.append(entropy)
        x.append(n_sphere_surfacearea(n = d - 2, r = T)/rho**((2-d)/d))
        
    plt.scatter(x, y_entropyList)
    plt.xlabel(r'$\frac{A}{\rho^{\frac{2-d}{d}}}$', fontsize = 20 )
    plt.ylabel(r'$s_{Boltz}$', fontsize = 20 )
    popt, pcov = curve_fit(linear, x, y_entropyList)
    xLinspace = np.linspace(min(x), max(x), 100)
    plt.plot(xLinspace, linear(xLinspace, *popt))
    plt.title(f's_Boltzmann in Dynamic in {d}d')
    plt.show()
    
    print(f'\n \n \n a_Boltzmann value for {d}d is {popt[0]}')
        
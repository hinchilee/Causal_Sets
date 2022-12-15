
import numpy as np
from causalsetfunctions import find_entropy, find_entropy_1molecule, convert_Harray_1molecule
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import json
import sys

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def linear(x, m, c): 
    return m*x + c 

d = {}
with open(path + 'H_i_Rindler.json') as f:
    d = json.load(f)
n_dict = {}
with open(path + 'H_i_n_Rindler.json') as f:
    n_dict = json.load(f)

rho_array = [int(ele) for ele in d.keys()] #[100, 500, 1000, 5000] 

d_normalised = dict() 
for rho in rho_array: 
    d_normalised[rho] = [ele/n_dict[str(rho)] for ele in d[str(rho)]]

# This script creates a plot and fits parameters for Rindler Horizon
# that checks that the scaling power is 1/2

#%%
# Obtaining entropy assuming counting the full range of moleules  

plt.figure(figsize = (12,8))
S = [find_entropy(d_normalised[rho]) for rho in rho_array]
plt.plot(rho_array, S, '.')

params, cov = curve_fit(linear, np.log10(rho_array), np.log10(S))
Nlinspace = np.linspace(min(rho_array), max(rho_array), 100)
plt.plot(Nlinspace, 10 ** linear(np.log10(Nlinspace), *params), label = 'Fit')
print(params)
print('err:', np.sqrt(np.diagonal(cov)))

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Sprinkling Density, $\rho$', fontsize = 20)
plt.ylabel('Entropy, s', fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.title('Counting molecules around Rindler Horizon (30 realisations)', fontsize = 20 )
plt.savefig('Plots/ Entropy_RindlerHorizon_Rel30_Molecules.png')
plt.show()

plt.clf()

#%%
# Obtaining entropy assuming counting 1-molecules only 

plt.figure(figsize = (12,8))
S = [find_entropy_1molecule(convert_Harray_1molecule(d_normalised[rho])) for rho in rho_array]
plt.plot(rho_array, S, '.')

params, cov = curve_fit(linear, np.log10(rho_array), np.log10(S))
Nlinspace = np.linspace(min(rho_array), max(rho_array), 100)
plt.plot(Nlinspace, 10 ** linear(np.log10(Nlinspace), *params), label = 'Fit')
print(params)
print('err:', np.sqrt(np.diagonal(cov)))

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Sprinkling Density, $\rho$', fontsize = 20)
plt.ylabel('Entropy, s', fontsize = 20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.title('Counting 1-molecules around Rindler Horizon (30 realisations)', fontsize = 20 )
plt.savefig('Plots/ Entropy_RindlerHorizon_Rel30_1Molecules.png')
plt.show()

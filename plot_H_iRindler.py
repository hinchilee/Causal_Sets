
import numpy as np
from causalsetfunctions import find_entropy, find_entropy_1molecule, convert_Harray_1molecule
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear(x, m, c): 
    return m*x + c 

rho_array = [100, 500, 1000, 5000]

df = pd.read_csv('H_iRindler.csv', index_col=0)

# This script creates a plot and fits parameters for Rindler Horizon
# that checks that the scaling power is 1/2

#%%
# Obtaining entropy assuming counting the full range of moleules  

plt.figure(figsize = (12,8))
S = [find_entropy(df[str(rho)].to_list()) for rho in rho_array]
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
S = [find_entropy_1molecule(convert_Harray_1molecule(df[str(rho)])) for rho in rho_array]
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

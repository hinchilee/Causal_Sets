import numpy as np
from causalsetfunctions import find_entropy
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear(x, m, c): 
    return m*x + c 

rho_array = [100, 500, 1000, 5000]

df = pd.read_csv('H_i.csv', index_col=0)
S = [find_entropy(df[str(rho)].to_list()) for rho in rho_array]
plt.plot(rho_array, S, '.')

params, cov = curve_fit(linear, np.log10(rho_array), np.log10(S))
Nlinspace = np.linspace(min(rho_array), max(rho_array), 100)
plt.plot(Nlinspace, 10 ** linear(np.log10(Nlinspace), *params), label = 'Fit')
print(params)

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\rho$')
plt.ylabel('Entropy')
plt.show()
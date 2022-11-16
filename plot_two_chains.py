import numpy as np
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

rho_array = [100, 500, 1000, 5000]
mean_array = []
sem_array = []

df = pd.read_csv('two_chains.csv')
for rho in rho_array:
    counts = df[str(rho)].tolist()
    mean_array.append(np.mean(counts))
    sem_array.append(stats.sem(counts))
    
def linear(x, m, c): 
    return m*x + c 

plt.errorbar(rho_array, mean_array, sem_array, fmt='o', label = 'data')
params, cov = curve_fit(linear, np.log10(rho_array), np.log10(mean_array))
Nlinspace = np.linspace(100,3e3)
plt.plot(Nlinspace, linear(Nlinspace,*params), label = 'fit')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel(r'Sprinkling Denisty, $\rho$')
plt.ylabel('Average Number of Links')
plt.show()
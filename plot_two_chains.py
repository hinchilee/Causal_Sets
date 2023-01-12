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

plt.errorbar(rho_array, mean_array, sem_array, fmt='o', label = 'Data')
params, cov = curve_fit(linear, np.log10(rho_array), np.log10(mean_array))
Nlinspace = np.linspace(min(rho_array), max(rho_array), 100)
plt.plot(Nlinspace, 10 ** linear(np.log10(Nlinspace), *params), label = 'Fit')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel(r'Sprinkling Denisty, $\rho$')
plt.ylabel('Average Number of Links')
plt.show()

plt.plot(rho_array, sem_array, 'o', label = 'Data')
params, cov = curve_fit(linear, np.log10(rho_array), np.log10(sem_array))
plt.plot(Nlinspace, 10 ** linear(np.log10(Nlinspace), *params), label = 'Fit')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel(r'Sprinkling Denisty, $\rho$')
plt.ylabel('Absolute Standard Error')
plt.show()

plt.plot(rho_array, np.divide(sem_array, mean_array), 'o', label = 'Data')
params, cov = curve_fit(linear, np.log10(rho_array), np.log10(np.divide(sem_array, mean_array)))
plt.plot(Nlinspace, 10 ** linear(np.log10(Nlinspace), *params), label = 'Fit')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel(r'Sprinkling Denisty, $\rho$')
plt.ylabel('Relative Standard Error')
plt.show()
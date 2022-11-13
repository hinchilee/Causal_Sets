import numpy as np
import scipy.stats as stats
import pandas as pd
import matplotlib.pyplot as plt

N_array = [100, 500, 1000, 5000]
mean_array = []
sem_array = []

df = pd.read_csv('two_chains.csv')
for N in N_array:
    counts = df[str(N)].tolist()
    mean_array.append(np.mean(counts))
    sem_array.append(stats.sem(counts))

plt.errorbar(N_array, mean_array, sem_array, fmt='o')
plt.xscale('log')
plt.xlabel('N')
plt.ylabel('Average Number of Links')
plt.show()
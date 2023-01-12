from causalset import CausalSet
import time
import numpy as np
import scipy.stats as stats
import pandas as pd

tic = time.time()

N_array = [100, 500, 1000, 5000, 10000, 50000]
mean_array = []
sem_array = []

for N in N_array:
    dimensions = [CausalSet(number_of_points=N, dimension=2).find_Myhreim_Meyer_dimension() for _ in range(5)]
    mean_array.append(np.mean(dimensions))
    sem_array.append(stats.sem(dimensions)[0])

df = pd.DataFrame({'N': N_array, 'mean': mean_array, 'sem': sem_array})
df.to_csv('dimension_estimates.csv', index=False)

toc = time.time()
print(f'Time elapsed is {toc - tic}')
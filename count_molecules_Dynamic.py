import numpy as np
from causalset import CausalSet
import time
import sys
import os
import csv
import pandas as pd

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def count_chains(N, T, Bounds):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    return CausalSet(sprinkling_density=N, dimension=4, BHtype = 'Dynamic', T=T, bounds = Bounds).find_molecules()

def main():
    tic = time.time()
    rho_array = [300]
    d_array = [4]
    T = 3
    for rho in rho_array:
        for dimension in d_array:
        # Number of realisations
            n = 1
            try:
                df = pd.read_csv(path + f'test_run_Dynamic_rho{rho}_{dimension}d.csv', names=['type', 'value'], header=None)
                min_time = max(df[df['type'] == 'min_time']['value'].min(), -1)
                min_distance = max(df[df['type'] == 'min_distance']['value'].min(), -0.5)
                max_distance = min(df[df['type'] == 'max_distance']['value'].max(), 1.5)

            except:
                min_time = -1
                print('No test run information!')
        
            boundsArray = [3 - min_distance, 3 + max_distance, T - min_time, T]
            for _i in range(n):
                print(f'\n realisation:{_i+1}, rho:{rho}, dimension:{dimension}')
                H = count_chains(rho, dimension, boundsArray)
                with open(path + f'H_Dynamic{dimension}d.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    writer.writerow([rho, H])
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()




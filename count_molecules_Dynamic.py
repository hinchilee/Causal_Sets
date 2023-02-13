import numpy as np
from causalset import CausalSet
import time
import multiprocessing as mp
import json
import sys
import os
import csv

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def count_chains(N, T, Bounds):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    return CausalSet(sprinkling_density=N, dimension=4, BHtype = 'Dynamic', T=T, bounds = Bounds).find_molecules()

def main():
    tic = time.time()
    with open(path + 'min_time.json') as f:
        min_times= json.load(f)
        
    rho_array = [300]
    d_array = [4]
    T = 3
    for rho in rho_array:
        for dimension in d_array:
        # Number of realisations
            n = 1
            try:
                min_time = min_times[f"{str(int(rho))}_{dimension}d"]
                if min_time == 10000: 
                    min_time = -1
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




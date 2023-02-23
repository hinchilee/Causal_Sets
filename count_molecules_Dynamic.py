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

def count_chains(N, d, Time, Bounds):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    return CausalSet(sprinkling_density=N, dimension=d, BHtype = 'Dynamic', T=Time, bounds = Bounds).find_molecules()

def main():
    tic = time.time()
    rho_array = [0.03]
    d_array = [4]
    T = 3
    for rho in rho_array:
        for dimension in d_array:
        # Number of realisations
            n = 100
            try:
                df = pd.read_csv(path + f'TestRun/test_run_Dynamic_rho{rho}_{dimension}d.csv', names=['type', 'value'], header=None)
                min_time = max(df[df['type'] == 'min_time']['value'].min()*1.1, -T)
                #min_time = -1.2
                #min_distance = 0
                #max_distance = 5.9
                min_distance = max(df[df['type'] == 'min_distance']['value'].min()*0.95, 0)
                max_distance = min(df[df['type'] == 'max_distance']['value'].max()*1.05, 2*T)

            except:
                raise ValueError ('No test run information!')
        
            #adding some leeway
            boundsArray = [min_distance, max_distance, T + min_time, T]
            for _i in range(n):
                print(f'\n realisation:{_i+1}, rho:{rho}, dimension:{dimension}')
                print('BoundsArray:', boundsArray)
                H = count_chains(rho, dimension, T, boundsArray)
                with open(path + f'H_Dynamic{dimension}d.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    writer.writerow([rho, H])
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()




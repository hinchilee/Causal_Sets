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

def count_chains(N, mintime, d = 4):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    boundsArray = np.array([[-0.5, 0.5] for i in range(d)])
    boundsArray[0][0] = mintime # no normalisation
    boundsArray[1][1] = 1.5
    #boundsArray[1][0] = mindistance 
    #boundsArray[1][1] = maxdistance
    
    return CausalSet(sprinkling_density=N, dimension=d, BHtype='Rindler', bounds = boundsArray).find_molecules()
  

def main():
    tic = time.time()
    rho_array = [10, 100]
    d_array = [2, 3]
    for rho in rho_array:
        for dimension in d_array:
        # Number of realisations
            n = 1
            try:
                df = pd.read_csv(path + f'test_run_Rindler_rho{rho}_{dimension}d.csv', names=['type', 'value'], header=None)
                min_time = max(df[df['type'] == 'min_time']['value'].min(), -1)
                min_distance = max(df[df['type'] == 'min_distance']['value'].min(), -0.5)
                max_distance = min(df[df['type'] == 'max_distance']['value'].max(), 1.5)

            except:
                min_time = -0.5
                print('No test run information!')
        
            for _i in range(n):
                print(f'\n realisation:{_i+1}, rho:{rho}, dimension:{dimension}')
                H = count_chains(rho, min_time, dimension)
                with open(path + f'H_Rindler{dimension}d.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    writer.writerow([rho, H])
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()

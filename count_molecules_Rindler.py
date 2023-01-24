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


def count_chains(N):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    return CausalSet(sprinkling_density=N, dimension=4, BHtype='Rindler').find_molecules()
  

def main():
    tic = time.time()

    rho_array = [1000, 3000]
    for rho in rho_array:
        # Number of realisations
        n = 2

        for _i in range(n):
            H = count_chains(rho)
            with open(path + 'H_Rindler.csv', 'a') as f:
                writer = csv.writer(f, lineterminator='\n')
                writer.writerow([rho, H])


if __name__ == "__main__":
    main()

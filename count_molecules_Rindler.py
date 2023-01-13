import numpy as np
from causalset import CausalSet
import time
import multiprocessing as mp
import json
import sys
import os

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def count_chains(N):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    return CausalSet(sprinkling_density=N, dimension=4, BHtype = 'Rindler').find_molecules()

def main():
    tic = time.time()

    rho_array = [1000] 
    with open(path + 'H_i_Rindler.json') as f:
        d = json.load(f)
    with open(path + 'H_i_n_Rindler.json') as f:
        n_dict = json.load(f)

    for rho in rho_array:
        # Number of realisations
        n = 3
        if str(rho) not in d.keys():
            d[str(rho)] = []
            n_dict[str(rho)] = 0
        n_dict[str(rho)] += n
        
        H_counts = [count_chains(rho) for _i in range(n)]

        max_H_i = max([len(H) for H in H_counts] + [len(d[str(rho)])])
        while len(d[str(rho)]) < max_H_i:
            d[str(rho)].append(0)
        for i in range(max_H_i):
            for H in H_counts:
                if len(H) > i:
                    d[str(rho)][i] += H[i]

        print(rho)
        toc = time.time()
        print(f'Time elapsed is {toc - tic}')
    
    print(d)
    print(n_dict)
    with open(path + 'H_i_Rindler.json', 'w', encoding='utf-8') as f:
        json.dump(d, f, ensure_ascii=False, indent=4)
    with open(path + 'H_i_n_Rindler.json', 'w', encoding='utf-8') as f:
        json.dump(n_dict, f, ensure_ascii=False, indent=4)

def reset(): 
    with open(path + 'H_i_Rindler.json', 'w', encoding='utf-8') as f:
        json.dump({}, f, ensure_ascii=False, indent=4)
    with open(path + 'H_i_n_Rindler.json', 'w', encoding='utf-8') as f:
        json.dump({}, f, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()

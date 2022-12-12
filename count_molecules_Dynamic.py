import numpy as np
from causalset import CausalSet
import time
import multiprocessing as mp
import json

def count_chains(N, T=5):
    return CausalSet(sprinkling_density=N, dimension=4, BHtype = 'Dynamic', T=T).find_molecules()

def main():
    tic = time.time()

    rho_array = [0.01, 0.03, 0.1, 0.3, 1] 
    d = {}
    with open('H_i_Dynamic.json') as f:
        d = json.load(f)
    n_dict = {}
    with open('H_i_n_Dynamic.json') as f:
        n_dict = json.load(f)

    for rho in rho_array:
        # Number of realisations
        n = 10
        if str(rho) not in d.keys():
            d[str(rho)] = []
            n_dict[str(rho)] = 0
        n_dict[str(rho)] += n

        pool = mp.Pool(mp.cpu_count() - 8)
        H_counts = pool.map(count_chains, [rho] * n)
        pool.close()

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
    with open('H_i_Dynamic.json', 'w', encoding='utf-8') as f:
        json.dump(d, f, ensure_ascii=False, indent=4)
    with open('H_i_n_Dynamic.json', 'w', encoding='utf-8') as f:
        json.dump(n_dict, f, ensure_ascii=False, indent=4)


def reset(): 
    with open('H_i_Dynamic.json', 'w', encoding='utf-8') as f:
        json.dump({}, f, ensure_ascii=False, indent=4)
    with open('H_i_n_Dynamic.json', 'w', encoding='utf-8') as f:
        json.dump({}, f, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()
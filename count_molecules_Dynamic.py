import numpy as np
from causalset import CausalSet
import time
import multiprocessing as mp

def count_chains(N, T=5):
    return CausalSet(sprinkling_density=N, dimension=4, BHtype = 'Dynamic', T=T).find_molecules()

def main():
    tic = time.time()

    rho_array = [0.01, 0.03, 0.1, 0.3, 1] 
    d = np.load('H_i_Dynamic.npy', allow_pickle=True).item()
    n_dict = np.load('H_i_n_Dynamic.npy', allow_pickle=True).item()

    for rho in rho_array:
        # Number of realisations
        n = 20
        if rho not in d.keys():
            d[rho] = []
            n_dict[rho] = 0
        n_dict[rho] += n

        pool = mp.Pool(mp.cpu_count() - 8)
        H_counts = pool.map(count_chains, [rho] * n)
        pool.close()

        max_H_i = max([len(H) for H in H_counts] + [len(d[rho])])
        while len(d[rho]) < max_H_i:
            d[rho].append(0)
        for i in range(max_H_i):
            for H in H_counts:
                if len(H) > i:
                    d[rho][i] += H[i]

        print(rho)
        toc = time.time()
        print(f'Time elapsed is {toc - tic}')
    
    print(d)
    print(n_dict)
    np.save('H_i_Dynamic.npy', d)
    np.save('H_i_n_Dynamic.npy', n_dict)

def reset(): 
    np.save('H_i_Dynamic.npy', {})
    np.save('H_i_n_Dynamic.npy', {})

if __name__ == "__main__":
    main()
from causalset import CausalSet
import time
import pandas as pd
import multiprocessing as mp

def count_chains(N):
    return CausalSet(sprinkling_density=N, dimension=4, DynamicBH=False).find_molecules()

def main():
    tic = time.time()

    rho_array = [100, 500, 1000, 5000]
    d = {}

    for rho in rho_array:
        # Number of realisations
        n = 50

        pool = mp.Pool(mp.cpu_count() - 8)
        H_counts = pool.map(count_chains, [rho] * n)
        pool.close()

        max_H_i = max([len(H) for H in H_counts])
        H_array = [0] * max_H_i
        for i in range(max_H_i):
            for H in H_counts:
                if len(H) > i:
                    H_array[i] += H[i] / n

        d[rho] = H_array

        print(rho)
        toc = time.time()
        print(f'Time elapsed is {toc - tic}')
    
    max_H_i = max([len(d[rho]) for rho in rho_array])
    for rho in rho_array:
        while len(d[rho]) < max_H_i:
            d[rho].append(0)

    df = pd.DataFrame(d)
    df.index += 1
    print(df)
    df.to_csv('H_i.csv')

if __name__ == "__main__":
    main()
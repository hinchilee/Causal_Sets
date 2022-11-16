from causalset import CausalSet
import time
import pandas as pd
import multiprocessing as mp

def count_chains(N):
    return CausalSet(number_of_points=N, dimension=4).find_molecules()

def main():
    tic = time.time()

    rho_array = [100, 500, 1000, 5000]#, 10000, 20000]
    d = {}

    for rho in rho_array:
        pool = mp.Pool(mp.cpu_count() - 4)
        # dimensions = [count_chains(N) for _ in range(5)]
        dimensions = pool.map(count_chains, [rho] * 20)
        pool.close() 
        d[rho] = dimensions

        print(rho)
        toc = time.time()
        print(f'Time elapsed is {toc - tic}')
    
    df = pd.DataFrame(d)
    df.to_csv('two_chains.csv')

if __name__ == "__main__":
    main()
from causalset import CausalSet
import time
import pandas as pd
import multiprocessing as mp

def count_chains(N):
    return CausalSet(number_of_points=N, dimension=4).find_molecules()

def main():
    tic = time.time()

    N_array = [100, 500, 1000, 5000]
    d = {}

    for N in N_array:
        pool = mp.Pool(mp.cpu_count() - 4)
        # dimensions = [count_chains(N) for _ in range(5)]
        dimensions = pool.map(count_chains, [N] * 10)
        pool.close() 
        d[N] = dimensions

        print(N)
        toc = time.time()
        print(f'Time elapsed is {toc - tic}')
    
    df = pd.DataFrame(d)
    df.to_csv('two_chains.csv')

if __name__ == "__main__":
    main()
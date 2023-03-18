import numpy as np
from causalset import CausalSet
from causalsetfunctions import generate_adjacency_matrix, count_subgraphs, check_connected_graph, compute_spacetimecuts_tube
import time
import sys
import os
import csv
import pandas as pd

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def count_chains(N, d, Time, moleculetype, Bounds):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    #lambda-molecules
    if moleculetype == 'lambda':
        return CausalSet(sprinkling_density=N, dimension=d, BHtype = 'Dynamic', sprinkling='Tube', T=Time, bounds = Bounds).find_molecules() 
    #v-moelcules
    elif moleculetype == 'v':
        c = CausalSet(sprinkling_density=N, dimension=d, BHtype = 'Dynamic', sprinkling='Tube', T=Time, bounds = Bounds)
        V = c.find_Vmolecules()
        try:
            Adjmatrix = generate_adjacency_matrix(c.VElementsLabelsList, c.CausalMatrix)
            SubGraphs = count_subgraphs(Adjmatrix)
            Connected = check_connected_graph(Adjmatrix)
            print(f'Subgraphs number: {SubGraphs}')
            print(f'Connected Graph: {Connected}')
            return V, SubGraphs, Connected
        except: 
            print('Subgraphs number: 0')
            print('Connected Graph: False')
            return V, 0, False
    
def main():
    tic = time.time()
    rho_array = [1000]
    d_array = [4]
    T = 1
    #moleculeType = 'v'
    moleculeType = 'v'
    #for rho in rho_array:
    for N_max in [10000]:
        for dimension in d_array:
        # Number of realisations
            n = 100
            # try:
            #     df = pd.read_csv(path + f'TestRun_T_{T}/test_run_Dynamic_rho{rho}_{dimension}d.csv', names=['type', 'value'], header=None)
            #     min_time = max(df[df['type'] == 'min_time']['value'].min()*1.1, -T)
            #     #min_time = -1.2
            #     #min_distance = 0
            #     #max_distance = 5.9
            #     min_distance = max(df[df['type'] == 'min_distance']['value'].min()*0.95, 0)
            #     max_distance = min(df[df['type'] == 'max_distance']['value'].max()*1.05, 2*T)

            # except:
            #     raise ValueError ('No test run information!')
        
            #adding some leeway
            #boundsArray = [min_distance, max_distance, T + min_time, T]
            #N_max = 10000
            boundsArray, rho = compute_spacetimecuts_tube(d = dimension, rho2 = 5000, N_max = N_max, b = 3)
                
            for _i in range(n):
                print(f'\n realisation:{_i+1}, rho:{rho}, dimension:{dimension}')
                print('BoundsArray:', boundsArray)
                print('N_max:', N_max)
                if moleculeType == 'lambda':
                    H = count_chains(rho, dimension, T, moleculeType, boundsArray)
                elif moleculeType == 'v':
                    H, Subgraphs, Connected = count_chains(rho, dimension, T, moleculeType, boundsArray)
                with open(path + f'H_Dynamic{dimension}d_{moleculeType}.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    if moleculeType == 'lambda':
                        writer.writerow([rho, H])
                    elif moleculeType == 'v':
                        writer.writerow([rho, H, Subgraphs, Connected])
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()




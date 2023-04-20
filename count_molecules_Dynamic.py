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
    if moleculetype == 'lambda' or moleculetype == 'v':
        c = CausalSet(sprinkling_density=N, dimension=d, BHtype = 'Dynamic', sprinkling='Tube', T=Time, bounds = Bounds)
        H, b =  c.find_molecules() 
        V, b2 = c.find_Vmolecules()
    #v-moelcules
    # elif moleculetype == 'v':
    #     c = CausalSet(sprinkling_density=N, dimension=d, BHtype = 'Dynamic', sprinkling='Tube', T=Time, bounds = Bounds)
    #     V, b2 = c.find_Vmolecules()
        try:
            Adjmatrix = generate_adjacency_matrix(c.VElementsLabelsList, c.CausalMatrix)
            SubGraphs = count_subgraphs(Adjmatrix)
            Connected = check_connected_graph(Adjmatrix)
            print(f'Subgraphs number: {SubGraphs}')
            print(f'Connected Graph: {Connected}')
            return H, b, V, SubGraphs, Connected, b2
        except: 
            print('Subgraphs number: 0')
            print('Connected Graph: False')
            return H, b, V, 0, False, b2
    
def main():
    tic = time.time()
    rho_array = [5000, 100000]
    d_array = [4]
    T = 1
    moleculeType = 'lambda'
    N_max = 23000
    #b = 4
    rho2 = 10
    #moleculeType = 'lambda'
    for dimension in d_array:
    #for N_max in [16000]:
        for b in [1.7, 1.5, 1.4, 1.7, 1.5, 1.4]:
        # Number of realisations
            n = 50
            
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
            boundsArray, rho = compute_spacetimecuts_tube(d = dimension, rho2 = rho2, N_max = N_max, b = b)
                
            for _i in range(n):
                print(f'\n realisation:{_i+1}, rho:{rho}, dimension:{dimension}')
                print('BoundsArray:', boundsArray)
                print('N_max:', N_max)
                if moleculeType == 'lambda' or moleculeType == 'v':
                    H , b, V, Subgraphs, Connected, b2 = count_chains(rho, dimension, T, moleculeType, boundsArray)
                with open(path + f'H_Dynamic{dimension}d_lambda.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    writer.writerow([rho, H, b])
                with open(path + f'H_Dynamic{dimension}d_v.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    writer.writerow([rho, V, Subgraphs, Connected, b2])
                
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()




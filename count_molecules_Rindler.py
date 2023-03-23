import numpy as np
from causalset import CausalSet
from causalsetfunctions import generate_adjacency_matrix, count_subgraphs, check_connected_graph
import time
import sys
import os
import csv
import pandas as pd
import json
from causalsetfunctions import compute_spacetimecuts_uniform_Rindler

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def count_chains(N, mintime, mindistance, maxdistance, moleculetype, boundsArray, d = 4):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    
    #lambda-molecules
    if moleculetype == 'lambda':
        return CausalSet(sprinkling_density=N, dimension=d, BHtype='Rindler', bounds = boundsArray).find_molecules()
    #v-moelcules
    elif moleculetype == 'v':
        c = CausalSet(sprinkling_density=N, dimension=d, BHtype='Rindler', bounds = boundsArray)
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
            print('Connected Graph: None')
            return V, 0, False
        
def main():
    tic = time.time()
    #rho_array = [1000]
    d_array = [4]
    moleculeType = 'lambda'
    #moleculeType = 'v'
    #for rho in rho_array:
    rho = 23000
    b = 1.9
    for N_max in [19000]:
        for dimension in d_array:
        # Number of realisations
            n = 50
            #rho = N_max 
            #if rho == 1000: 
             #   rho -= 100
            # try:
            #     df = pd.read_csv(path + f'TestRun/test_run_Rindler_rho{rho}_{dimension}d.csv', names=['type', 'value'], header=None)
            #     with open(path + 'min_time.json') as f:
            #         dfTime = json.load(f)
            #     min_time = max(dfTime[f"{rho}_{dimension}d"], -1)
            #     min_distance = max(df[df['type'] == 'min_distance']['value'].min()*0.95, 0)
            #     max_distance = min(df[df['type'] == 'max_distance']['value'].max()*1.05, 2)

            # except:
                
            #     raise ValueError ('No test run information!')
            
            boundsArray, adjusted_rho, l, adjusted_l = compute_spacetimecuts_uniform_Rindler(d = dimension, rho = rho, N_max = N_max, b= b)
        
            for _i in range(n):
                print(f'\n realisation:{_i+1}, rho:{adjusted_rho}, dimension:{dimension}')
                print('BoundsArray:\n', boundsArray)
                print('N_max:', N_max)
                print(f'l: {l}, adjustedl: {adjusted_l}')
                if moleculeType == 'lambda':
                    H = count_chains(adjusted_rho, 0, 0, 0, moleculeType, boundsArray, dimension)
                elif moleculeType == 'v':
                    H, Subgraphs, Connected = count_chains(adjusted_rho, 0, 0, 0, moleculeType, boundsArray, dimension)
                with open(path + f'H_Rindler{dimension}d_{moleculeType}.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    if moleculeType == 'lambda':
                        writer.writerow([adjusted_rho, H])
                    elif moleculeType == 'v':
                        writer.writerow([adjusted_rho, H, Subgraphs, Connected])
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()

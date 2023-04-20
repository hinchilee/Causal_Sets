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
    if moleculetype == 'lambda' or moleculetype == 'v':
        c = CausalSet(sprinkling_density=N, dimension=d, BHtype='Rindler', bounds = boundsArray)
        H, b = c.find_molecules()
    #v-moelcule
        V, b2 = c.find_Vmolecules()
        try:
            Adjmatrix = generate_adjacency_matrix(c.VElementsLabelsList, c.CausalMatrix)
            SubGraphs = count_subgraphs(Adjmatrix)
            Connected = check_connected_graph(Adjmatrix)
            print(f'Subgraphs number: {SubGraphs}')
            print(f'Connected Graph: {Connected}')
            return H, b, V, SubGraphs, Connected, b2
        except: 
            print('Subgraphs number: 0')
            print('Connected Graph: None')
            return H, b, V, 0, False, b2
        
def main():
    tic = time.time()
    #d_array = [2]
    d_array = [4]
    moleculeType = 'v'
    #N_max = 1000
    N_max = 7000
    #for N_max in [19000]:
    for dimension in d_array:
        if dimension == 4: 
            rho_array = [1000, 30000]
        for rho in rho_array:
        # Number of realisations
            n = 100
            if rho == 30000: 
                b = 2.8
                N_max = 10000
            elif rho == 1000:
                b = 3
                N_max = 8000
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
                if moleculeType == 'lambda' or moleculeType == 'v':
                    H , b, V, Subgraphs, Connected, b2 = count_chains(adjusted_rho, 0, 0, 0, moleculeType, boundsArray, dimension)
                with open(path + f'H_Rindler{dimension}d_lambda.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    writer.writerow([adjusted_rho, H, b])
                with open(path + f'H_Rindler{dimension}d_v.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    writer.writerow([adjusted_rho, V, Subgraphs, Connected, b2])
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()

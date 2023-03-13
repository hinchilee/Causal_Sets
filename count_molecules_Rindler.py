import numpy as np
from causalset import CausalSet
from causalsetfunctions import generate_adjacency_matrix, count_subgraphs, check_connected_graph
import time
import sys
import os
import csv
import pandas as pd
import json

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def count_chains(N, mintime, mindistance, maxdistance, moleculetype, d = 4):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    boundsArray = np.array([[-0.5, 0.5] for i in range(d)])
    boundsArray[0][0] = 0.5 + mintime 
    boundsArray[1][1] = 1.5
    boundsArray[1][0] = mindistance - 0.5
    boundsArray[1][1] = maxdistance - 0.5
    
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
    rho_array = [100, 300, 1000, 3000, 10000]
    d_array = [2,3,4]
    moleculeType = 'v'
    #moleculeType = 'lambda'
    for rho in rho_array:
        for dimension in d_array:
        # Number of realisations
            n = 100
            try:
                df = pd.read_csv(path + f'TestRun/test_run_Rindler_rho{rho}_{dimension}d.csv', names=['type', 'value'], header=None)
                with open(path + 'min_time.json') as f:
                    dfTime = json.load(f)
                min_time = max(dfTime[f"{rho}_{dimension}d"], -1)
                min_distance = max(df[df['type'] == 'min_distance']['value'].min()*0.95, 0)
                max_distance = min(df[df['type'] == 'max_distance']['value'].max()*1.05, 2)

            except:
                
                raise ValueError ('No test run information!')
        
            for _i in range(n):
                print(f'\n realisation:{_i+1}, rho:{rho}, dimension:{dimension}')
                if moleculeType == 'lambda':
                    H = count_chains(rho, min_time, min_distance, max_distance, moleculeType, dimension)
                elif moleculeType == 'v':
                    H, Subgraphs, Connected = count_chains(rho, min_time, min_distance, max_distance, moleculeType, dimension)
                with open(path + f'H_Rindler{dimension}d_{moleculeType}.csv', 'a') as f:
                    writer = csv.writer(f, lineterminator='\n')
                    if moleculeType == 'lambda':
                        writer.writerow([rho, H])
                    elif moleculeType == 'v':
                        writer.writerow([rho, H, Subgraphs, Connected])
                
    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()

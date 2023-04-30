import numpy as np
from causalset import CausalSet
from causalsetfunctions import generate_adjacency_matrix, count_subgraphs, check_connected_graph
import time
import sys
import os
import csv
import pandas as pd
import json
from causalsetfunctions import compute_spacetimecuts_uniform_Rindler, n_sphere_surfacearea

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'


def count_chains(N, mintime, mindistance, maxdistance, moleculetype, boundsArray, d=4):
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))

    # lambda-molecules
    if moleculetype == 'lambda' or moleculetype == 'v':
        c = CausalSet(sprinkling_density=N, dimension=d,
                      BHtype='Rindler', bounds=boundsArray)
        H, b = c.find_molecules()
    # v-moelcule
        V, b2 = c.find_Vmolecules()
        try:
            Adjmatrix = generate_adjacency_matrix(
                c.VElementsLabelsList, c.CausalMatrix)
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
    d_array = [3, 4]
    moleculeType = 'lambda'
    N_max = 10000
    c = 4
    rho2 = 10
    # for N_max in [19000]:
    for dimension in d_array:
        if dimension == 3:

            for N_max in [9000, 7000]:
                # Number of realisations
                n = 100
                c = 4

                boundsArray, adjusted_rho, l = compute_spacetimecuts_uniform_Rindler(
                    d=dimension, rho2=rho2, N_max=N_max, b=c)

                for _i in range(n):
                    print(
                        f'\n realisation:{_i+1}, rho:{adjusted_rho}, dimension:{dimension}')
                    print('BoundsArray:\n', boundsArray)
                    print('N_max:', N_max)
                    print(f'l: {l}')
                    if moleculeType == 'lambda' or moleculeType == 'v':
                        H, b, V, Subgraphs, Connected, b2 = count_chains(
                            adjusted_rho, 0, 0, 0, moleculeType, boundsArray, dimension)
                    with open(path + f'H_Rindler{dimension}d_lambda_additional.csv', 'a') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow([adjusted_rho, H, b])
                    with open(path + f'H_Rindler{dimension}d_v_additional.csv', 'a') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow(
                            [adjusted_rho, V, Subgraphs, Connected, b2])

        elif dimension == 4:

            for c in [3.8, 4.1]:
                # Number of realisations
                n = 100
                N_max = 10000

                boundsArray, adjusted_rho, l = compute_spacetimecuts_uniform_Rindler(
                    d=dimension, rho2=rho2, N_max=N_max, b=c)

                for _i in range(n):
                    print(
                        f'\n realisation:{_i+1}, rho:{adjusted_rho}, dimension:{dimension}')
                    print('BoundsArray:\n', boundsArray)
                    print('N_max:', N_max)
                    print(f'l: {l}')
                    if moleculeType == 'lambda' or moleculeType == 'v':
                        H, b, V, Subgraphs, Connected, b2 = count_chains(
                            adjusted_rho, 0, 0, 0, moleculeType, boundsArray, dimension)
                    with open(path + f'H_Rindler{dimension}d_lambda_additional.csv', 'a') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow([adjusted_rho, H, b])
                    with open(path + f'H_Rindler{dimension}d_v_additional.csv', 'a') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow(
                            [adjusted_rho, V, Subgraphs, Connected, b2])

    toc = time.time()
    print(f'Time taken is {toc - tic}')


if __name__ == "__main__":
    main()

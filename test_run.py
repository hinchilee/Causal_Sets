from causalset import CausalSet
import sys
import numpy as np 
import os
import csv

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def main():
    for BHtype in ['Rindler', 'Dynamic']:
        if BHtype == 'Rindler':
            rho_array = [10, 100]
        elif BHtype == 'Dynamic':
            rho_array = [1, 1.5]
        for rho in rho_array:
            for dimension in [2, 3]:
                print('sprinkling density', rho)
                print('dimension:', dimension)
                for _i in range(1):
                    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
                    boundsArray = np.array([[-0.5, 0.5] for i in range(dimension)])
                    boundsArray[1][1] = 1.5
                    c = CausalSet(sprinkling_density=rho, dimension=dimension, BHtype = BHtype, bounds = boundsArray)
                    c.find_molecules()

                    with open(path + f'test_run_{BHtype}_rho{rho}_{dimension}d.csv', 'a') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow(['min_time', c.min_time])
                        writer.writerow(['min_distance', c.min_distance])
                        writer.writerow(['max_distance', c.max_distance])

if __name__ == "__main__":
    main()
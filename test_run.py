from causalset import CausalSet
import sys
import numpy as np 
import os
import csv
import time
import json

path = 'TestRun/'
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def main():
    tic = time.time()
    for BHtype in ['Rindler']:
        if BHtype == 'Rindler':
            rho_array = [100, 300, 1000, 3000, 10000]
        elif BHtype == 'Dynamic':
            rho_array = [0.1, 0.3, 1, 3] 
            T = 3 
        for rho in rho_array:
            for dimension in [2]:
                #iterations
                for _i in range(40):
                    print('sprinkling density', rho)
                    print('dimension:', dimension)
                    print('realisation', _i)
                    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
                    if BHtype == 'Rindler':
                        path2 = ''
                        with open(path2 + 'min_time.json') as f:
                            dfTime = json.load(f)
                        min_time = max(dfTime[f"{rho}_{dimension}d"], -1)
                        boundsArray = np.array([[-0.5, 0.5] for i in range(dimension)])
                        boundsArray[1][1] = 1.5
                        boundsArray[0][0] = 0.5 + min_time 
                        c = CausalSet(sprinkling_density=rho, dimension=dimension, BHtype = BHtype, bounds = boundsArray)
                        c.find_molecules()
                    elif BHtype == 'Dynamic':
                        boundsArray = [0, 2*T, 0, T]
                        c = CausalSet(sprinkling_density=rho, dimension=dimension, BHtype = 'Dynamic', T=T, bounds = boundsArray)
                        c.find_molecules()
                    
                    with open(path + f'test_run_{BHtype}_rho{rho}_{dimension}d.csv', 'a') as f:
                        writer = csv.writer(f, lineterminator='\n')
                        writer.writerow(['min_time', c.min_time])
                        writer.writerow(['min_distance', c.min_distance])
                        writer.writerow(['max_distance', c.max_distance])
    toc= time.time() 
    print(f'Time taken is {toc - tic}')
if __name__ == "__main__":
    main()
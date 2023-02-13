from causalset import CausalSet
import sys
import json
import numpy as np 

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def main():
    with open(path + 'min_timeRindler.json') as f:
        d = json.load(f)

    for rho in [3000, 10000]:
        for dimension in [2,3]:
            print('sprinkling density', rho)
            print('dimension:', dimension)
            if (f"{str(rho)}_{dimension}d") not in d.keys():
                d[f"{str(rho)}_{dimension}d"] = []
        
            for i in range(5):
                boundsArray = np.array([[-0.5, 1.5] for i in range(dimension)])
                boundsArray[0][1] = 0.5
                c = CausalSet(sprinkling_density=rho, dimension=dimension, BHtype = 'Rindler', bounds = boundsArray)
                c.find_molecules()
                d[f"{str(rho)}_{dimension}d"].append(c.min_time)
            
            with open(path + 'min_timeRindler.json', 'w', encoding='utf-8') as f:
                json.dump(d, f, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()
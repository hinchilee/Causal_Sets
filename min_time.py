from causalset import CausalSet
import sys
import json

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def main():
    with open(path + 'min_time.json') as f:
        d = json.load(f)

    for rho in [1000, 3000]:
        for dimension in [2,3,4]:
            print('sprinkling density', rho)
            print('dimension:', dimension)
            if (f"{str(rho)}_{dimension}d") not in d.keys():
                d[f"{str(rho)}_{dimension}d"] = []
        
            for i in range(5):
                c = CausalSet(sprinkling_density=rho, dimension=dimension, BHtype = 'Rindler')
                c.find_molecules()
                d[f"{str(rho)}_{dimension}d"].append(c.min_time)
            
            with open(path + 'min_time.json', 'w', encoding='utf-8') as f:
                json.dump(d, f, ensure_ascii=False, indent=4)

if __name__ == "__main__":
    main()
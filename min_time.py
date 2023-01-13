from causalset import CausalSet
import sys
import json

path = ''
if len(sys.argv) > 1:
    path = sys.argv[1] + '/'

def main():
    with open(path + 'min_time.json') as f:
        d = json.load(f)

    for density in [1000]:
        print(density)
        if str(density) not in d.keys():
            d[str(density)] = []

        for i in range(5):
            c = CausalSet(sprinkling_density=density, dimension=4, BHtype = 'Rindler')
            c.find_molecules()
            d[str(density)].append(c.min_time)
        
        with open(path + 'min_time.json', 'w', encoding='utf-8') as f:
            json.dump(d, f, ensure_ascii=False, indent=4)


if __name__ == "__main__":
    main()
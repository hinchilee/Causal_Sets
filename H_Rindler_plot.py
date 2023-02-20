import pandas as pd
import matplotlib.pyplot as plt

def total_links(H):
    return sum([(i + 1) * H[i] for i in range(len(H))])

d_array = [2,3,4]
for d in d_array: 
    df = pd.read_csv(f'H_Rindler{d}d.csv', names=['rho', 'H'], header=None)
    print(df)
    rho_array = df['rho'].unique()
    rho_array.sort()
    for rho in rho_array:
        a_array = []
        H_array = [[int(x) for x in a.replace('[', '').replace(']', '').split(', ') if x != ''] for a in df[df['rho'] == rho]['H'].values]
        for H in H_array:
            a_array.append(rho**((2-d)/d)*(total_links(H)))

        plt.hist(a_array, density=True, label=rho, histtype='step', range=(0, 0.5))
    plt.legend(title=r'$\rho$')
    plt.title(f'({d-1} + 1) dimensions')
    plt.show()
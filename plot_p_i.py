import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
plt.rc('font', family='Arial', size=20)

bar_index = np.arange(6)
width = 0.35
d_array = [2, 3, 4]
for BHtype in ['Rindler', 'Dynamic']:
    df = pd.concat([pd.read_csv(f'H_{BHtype}{d}d.csv', names=['rho', 'H'], header=None) for d in d_array])
    H_dict = {}
    N = df.shape[0]
    for index, row in df.iterrows():
        H = [int(x) if x != '' else 0 for x in row['H'].replace('[', '').replace(']', '').split(', ')]
        for i in range(len(H)):
            if (i + 1) in H_dict:
                H_dict[i + 1].append(H[i])
            else:
                H_dict[i + 1] = [H[i]]

    mean_array = []
    std_array = []
    for H_i in H_dict.values():
        H_i = np.pad(H_i, (0, N - len(H_i)), 'constant')
        mean_array.append(np.mean(H_i))
        std_array.append(stats.sem(H_i))
    std_array /= np.sum(mean_array)
    mean_array /= np.sum(mean_array)
    if BHtype == 'Rindler':
        padding = 0
    else:
        padding = width
    plt.bar(bar_index + padding, mean_array, width, yerr=std_array, label=BHtype, capsize=10)
plt.ylabel(r'$p_i$')
plt.xlabel(r'$i$-molecule')
plt.xticks(bar_index + width / 2, ['1', '2', '3', '4', '5', '6'])
plt.legend()
plt.tight_layout()
plt.savefig(f'Plots/p_i.png', dpi=300)
# plt.show()
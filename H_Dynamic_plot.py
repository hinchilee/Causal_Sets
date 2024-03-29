import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from causalsetfunctions import n_sphere_surfacearea

plt.rc('font', family='Arial', size=20)

# theoretical_a = {
#     2: 1/3,
#     3: 0.218853,
#     4: 0.173205
# }

def H(s):
    return [int(x) for x in s.replace('[', '').replace(']', '').split(', ') if x != '']

def total_links(H):
    return sum([(i + 1) * H[i] for i in range(len(H))])

d_array = [2, 3, 4]
# d_array = [2]
for d in d_array: 
    df = pd.read_csv(f'H_Dynamic{d}d.csv', names=['rho', 'H'], header=None)
    df['a'] = df.apply(lambda row: row['rho']**((2-d)/d)*(total_links(H(row['H']))), axis=1) / n_sphere_surfacearea(n = d - 2, r = 3)
    sns.catplot(x='rho', y='a', data=df, height=5, aspect=1.5, kind='violin', cut=0)
    # plt.axhline(theoretical_a[d], color='r', linestyle='--')
    plt.xlabel(r'$\rho$')
    plt.ylabel('a')
    plt.title(f'Dynamic ({d-1} + 1) dimensions')
    plt.tight_layout()
    plt.savefig(f'Plots/H_Dynamic{d}d.png', dpi=300)
    # plt.show()
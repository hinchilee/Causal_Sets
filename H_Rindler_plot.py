import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def H(s):
    return [int(x) for x in s.replace('[', '').replace(']', '').split(', ') if x != '']

def total_links(H):
    return sum([(i + 1) * H[i] for i in range(len(H))])

d_array = [2, 3, 4]
for d in d_array: 
    df = pd.read_csv(f'H_Rindler{d}d.csv', names=['rho', 'H'], header=None)
    df['a'] = df.apply(lambda row: row['rho']**((2-d)/d)*(total_links(H(row['H']))), axis=1)
    sns.catplot(x='rho', y='a', data=df, height=5, aspect=2, kind='violin')
    plt.xlabel(r'$\rho$')
    plt.ylabel('a')
    plt.title(f'({d-1} + 1) dimensions')
    plt.show()
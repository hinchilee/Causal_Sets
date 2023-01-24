import pandas as pd
import itertools

df = pd.read_csv('H_Rindler.csv', names=['rho', 'H'], header=None)
print(df)
rho_array = df['rho'].unique()
for rho in rho_array:
    print(rho)
    print(df[df['rho'] == rho].shape[0])
    print([sum(x) for x in itertools.zip_longest(*[[int(x) for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)])
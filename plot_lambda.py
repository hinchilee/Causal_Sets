# import matplotlib
# matplotlib.use('pgf')
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })

import matplotlib.pyplot as plt

import numpy as np
import scipy.stats as stats
import pandas as pd
import itertools
from causalsetfunctions import n_sphere_surfacearea
import matplotlib

matplotlib.rcParams.update(matplotlib.rcParamsDefault)
params = {'font.size': 25,
          'font.family': 'Times New Roman',
          'axes.labelsize': 29,
          'legend.fontsize': 23,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'figure.figsize': [8.5, 6.5],
          'axes.prop_cycle': plt.cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color']
                                        + ['magenta'])
          }
plt.rcParams.update(params)


def count_links(a):
    links = 0
    for i in range(len(a)):
        links += (i + 1) * int(a[i])

    return links


def entropy(a):
    H = np.array([int(i) for i in a])
    N = np.sum(H)
    p = H / N
    S = N * np.sum([p_i * np.log(1 / p_i) for p_i in p if p_i != 0])

    return S


def flatten(l):
    return list(itertools.chain.from_iterable(l))


horizon_types = ['Rindler', 'Dynamic']
d_array = [2, 3, 4]
for horizon in horizon_types:
    print(horizon)
    for d in d_array:
        df = pd.read_csv(f'H_{horizon}{d}d_lambda_clean.csv', names=[
                         'rho', 'H', 'b'], header=None)

        rho_array = df['rho'].unique()
        rho_array.sort()

        H_links = []
        H_links_err = []
        A = []
        links_i = []

        for rho in rho_array:
            links = [entropy(a.replace('[', '').replace(']', '').split(
                ', ')) if a != '[]' else 0 for a in df[df['rho'] == rho]['H'].values]
            links_i.append(links)
            H_links.append(np.mean(links))
            H_links_err.append(stats.sem(links))

            l = rho**(-1/d)

            if horizon == 'Rindler':
                A.append(1 / (l ** (d - 2)))
            elif horizon == 'Dynamic':
                A.append(n_sphere_surfacearea(n=d-2, r=1) / (l ** (d - 2)))

        if d == 2:
            links = flatten(links_i)
            if horizon == 'Rindler':
                a = np.mean(links)
                a_err = stats.sem(links)
            elif horizon == 'Dynamic':
                a = np.mean(links) / 2
                a_err = stats.sem(links) / 2
            print(f'{round(a, 3)} +- {round(a_err, 3)}')

        else:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                A, H_links)
            print(f'{round(slope, 3)} +- {round(std_err, 3)}')

            plt.errorbar(A, H_links, fmt='o', yerr=H_links_err,
                         label='Data', capsize=5, ecolor='black')
            plt.plot(A, [slope * a + intercept for a in A],
                     label='Linear Fit', zorder=0, color='red')
            plt.ylabel(r'$S_{Boltz}$')
            if d == 3:
                plt.xlabel(r'$A/{\ell^{d-2}}$')
            elif d == 4:
                plt.xlabel(r'$A/{\ell^{d-2}}$')
            plt.legend()
            plt.tight_layout()
            plt.savefig(
                fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\Lambda_{horizon}_{d}d.png', dpi=300, bbox_inches='tight', transparent=True)

            plt.show()
            # plt.savefig(f'link_{d}d_{horizon}.pgf')
            # plt.clf()

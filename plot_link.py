# import matplotlib
# matplotlib.use('pgf')
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })

import numpy as np
import matplotlib
from scipy.optimize import curve_fit
from causalsetfunctions import n_sphere_surfacearea
import itertools
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

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


def flatten(l):
    return list(itertools.chain.from_iterable(l))


horizon_types = ['Rindler', 'Dynamic']
d_array = [2, 3, 4]
# horizon_types = ['Rindler', 'Dynamic']
# d_array = [2, 3, 4]
for horizon in horizon_types:
    print(horizon)
    for d in d_array:
        df = pd.read_csv(f'H_{horizon}{d}d_lambda_clean.csv', names=[
                         'rho', 'H', 'b'], header=None)

        rho_array = df['rho'].unique()
        rho_array.sort()

        H_links = []
        H_links_err = []
        links_i = []
        l = []

        for rho in rho_array:
            links = [count_links(a.replace('[', '').replace(']', '').split(
                ', ')) if a != '[]' else 0 for a in df[df['rho'] == rho]['H'].values]
            links_i.append(links)
            H_links.append(np.mean(links))
            H_links_err.append(stats.sem(links))

            l.append(rho**(-1/d))

        def area(l, k):
            if horizon == 'Rindler':
                return [1 for _ in l]
            elif horizon == 'Dynamic':
                return [n_sphere_surfacearea(n=d-2, r=1) * np.power(1 - (k * l_i), d - 2) for l_i in l]

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
            if horizon == 'Rindler':
                A = area(l, 0) / np.power(l, d-2)
                slope, intercept, r_value, p_value, std_err = stats.linregress(
                    A, H_links)
                print(f'{round(slope, 3)} +- {round(std_err, 3)}')
                plt.errorbar(A, H_links, fmt='o', yerr=H_links_err,
                             label='Data', capsize=5, ecolor='black')
                plt.plot(A, [slope * a + intercept for a in A],
                         label='Linear Fit', zorder=0, color='red')

            elif horizon == 'Dynamic':
                l = np.array(l)

                def H(l, a, k):
                    if d == 4:
                        b_3 = -0.0209261
                    else:
                        b_3 = 0
                    return ((a + ((b_3/np.sqrt(2)) * (d - 2) * l)) * area(l, k)) / np.power(l, d-2)
                    # return (a * area(l, k)) / np.power(l, d-2)
                if d == 3:
                    p0 = [0.219, 1]
                elif d == 4:
                    p0 = [0.173, 1]
                popt, pcov = curve_fit(
                    H, l, H_links, p0=p0, sigma=H_links_err, absolute_sigma=True)
                print(popt)
                print(np.sqrt(np.diag(pcov)))

                A = area(l, popt[1]) / np.power(l, d-2)
                plt.errorbar(A, H_links, fmt='o', yerr=H_links_err,
                             label='Data', capsize=5, ecolor='black')
                plt.plot(A, H(l, *popt), label='Linear Fit',  color='red')

            plt.ylabel(r'$\langle H_{links} \rangle$')
            if d == 3:
                plt.xlabel(r'$A/{\ell^{d-2}}$')
            elif d == 4:
                plt.xlabel(r'$A/{\ell^{d-2}}$')
            plt.legend()
            plt.tight_layout()
            plt.savefig(
                fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\Links_{horizon}_{d}d.png', dpi=300, bbox_inches='tight', transparent=True)
            plt.show()
            # plt.savefig(f'link_{d}d_{horizon}.pgf')
            # plt.clf()
            # plt.clf()
            # plt.clf()

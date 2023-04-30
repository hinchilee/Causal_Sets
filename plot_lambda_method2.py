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
from scipy.optimize import curve_fit

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


def entropy(N, p_iDict):

    S = N * np.sum([p_i * np.log(1 / p_i)
                   for p_i in p_iDict.values() if p_i != 0])

    return S


def flatten(l):
    return list(itertools.chain.from_iterable(l))


horizon_types = ['Rindler', 'Dynamic']
d_array = [3, 4]
for horizon in horizon_types:
    print(horizon)
    for d in d_array:
        df = pd.read_csv(f'H_{horizon}{d}d_lambda_clean.csv', names=[
                         'rho', 'H', 'b'], header=None)

        rho_array = df['rho'].unique()
        rho_array.sort()

        H = []
        N = []
        p_i = {}
        p_i_err = {}
        p_i_mean = {}
        A = []
        s_err_rho = []
        s_rho = []

        for rho in rho_array:
            H = [a.replace('[', '').replace(']', '').split(
                ', ') if a != '[]' else 0 for a in df[df['rho'] == rho]['H'].values]

            def find_max_length_list(list):
                list_len = [len(i) for i in list]
                return(max(list_len))

            max_lambda = find_max_length_list(H)

            # Initialise p_i {}
            for i in range(1, max_lambda + 1):
                p_i[i] = []

            for listH in H:
                # Extend all lists inside H to length of max_lambda
                while len(listH) < max_lambda:
                    listH.append('0')
                arrayH = np.array([int(ele) for ele in listH])
                N.append(sum([ele for ele in arrayH]))

                for i, H in enumerate(arrayH):

                    p_i[i+1].append(H/np.sum(arrayH))

            for i in p_i.keys():
                p_i_err[i] = stats.sem(p_i[i])
                p_i_mean[i] = np.mean(p_i[i])

            N_mean = np.mean(N)
            # print(N_mean)
            N_err = stats.sem(N)

            def error_propagate(p_i_mean, p_i_err, N_mean, N_err):
                s_err2 = 0
                # print(p_i_mean)
                for i in p_i_mean.keys():
                    if p_i_mean[i] != 0:  # Avoid log(0) error
                        s_err2 += ((np.log(p_i_mean[i]) + 1)
                                   )**2 * p_i_err[i]**2
                        s_err2 += (p_i_mean[i] *
                                   np.log(p_i_mean[i]))**2 * N_err**2
                    else:
                        pass

                # print(s_err2)

                return s_err2**0.5

            s_err = error_propagate(p_i_mean, p_i_err, N_mean, N_err)
            s_err_rho.append(s_err)
            s = entropy(N_mean, p_i_mean)
            s_rho.append(s)

            l = rho**(-1/d)

            if horizon == 'Rindler':
                A.append(1 / (l ** (d - 2)))
            elif horizon == 'Dynamic':
                A.append(n_sphere_surfacearea(n=d-2, r=1) / (l ** (d - 2)))

        if d == 2:
            # links = flatten(links_i)
            # if horizon == 'Rindler':
            #     a = np.mean(links)
            #     a_err = stats.sem(links)
            # elif horizon == 'Dynamic':
            #     a = np.mean(links) / 2
            #     a_err = stats.sem(links) / 2
            # print(f'{round(a, 3)} +- {round(a_err, 3)}')
            pass

        else:
            def linear(x, m, c):
                return m*x + c

            popt, pcov = curve_fit(linear,
                                   A, s_rho, sigma=s_err_rho, absolute_sigma=True)
            print(f'{round(popt[0], 3)} +- {round(np.sqrt(pcov[0,0]), 3)}')

            plt.errorbar(A, s_rho, fmt='o', yerr=s_err_rho,
                         label='Data', capsize=5, ecolor='black')
            plt.plot(A, linear(np.array(A), *popt),
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

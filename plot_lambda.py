# import matplotlib
# matplotlib.use('pgf')
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })

import matplotlib.pyplot as plt

def set_size(width_pt, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to sit nicely in our document.

    Parameters
    ----------
    width_pt: float
            Document width in points
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = 3 / 4

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)

# plt.figure(figsize=set_size(483.69687 * 0.45))

import numpy as np
import scipy.stats as stats
import pandas as pd
import itertools
from causalsetfunctions import n_sphere_surfacearea

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

# horizon_types = ['Rindler', 'Dynamic']
# d_array = [2]
horizon_types = ['Rindler', 'Dynamic']
d_array = [2, 3, 4]

for horizon in horizon_types:
    print(horizon)
    for d in d_array:
        df = pd.read_csv(f'H_{horizon}{d}d_lambda_clean.csv', names=['rho', 'H', 'b'], header=None)

        rho_array = df['rho'].unique()
        rho_array.sort()

        H_links = []
        H_links_err = []
        A = []
        links_i = []

        H_i = {}
        N_array = []
        if d == 2:
            for rho in rho_array:
                H_array = [a.replace('[', '').replace(']', '').split(', ') for a in df[df['rho'] == rho]['H'].values if a != '[]']
                for H in H_array:
                    H = [int(i) for i in H]
                    N_array.append(sum(H))
                    for i in range(len(H)):
                        if H[i] != 0:
                            if (i + 1) in H_i:
                                H_i[i + 1].append(H[i])
                            else:
                                H_i[i + 1] = [H[i]]
            for i in H_i:
                H_i[i] += [0] * (len(H_array) - len(H_i[i]))
            n_i = np.array([sum(H_i[i]) for i in H_i])
            n_error_i = np.array([stats.sem(H_i[i]) for i in H_i])
            N = sum(n_i)
            # print(n_i)
            p_i = n_i / N
            p_error_i = n_error_i
            # print(p_i)
            S = N/len(df) * np.sum([p * np.log(1 / p) for p in p_i if p != 0])
            N_error = stats.sem(N_array)
            S_error = np.sqrt(np.sum([((np.log(p_i[i]) + 1) ** 2) * (p_error_i[i] ** 2) + ((p_i[i] * np.log(p_i[i])) ** 2) * (N_error ** 2) for i in range(len(p_i))]))

            if horizon == 'Rindler':
                print(f'{round(S, 3)} +- {round(S_error, 3)}')
            elif horizon == 'Dynamic':
                print(f'{round(S / 2, 3)} +- {round(S_error / 2, 3)}')

            continue


        for rho in rho_array:
            H_i = {}
            H_array = [a.replace('[', '').replace(']', '').split(', ') for a in df[df['rho'] == rho]['H'].values if a != '[]']
            N_array = []
            for H in H_array:
                H = [int(i) for i in H]
                N_array.append(sum(H))
                for i in range(len(H)):
                    if H[i] != 0:
                        if (i + 1) in H_i:
                            H_i[i + 1].append(H[i])
                        else:
                            H_i[i + 1] = [H[i]]
            for i in H_i:
                H_i[i] += [0] * (len(H_array) - len(H_i[i]))
            n_i = np.array([sum(H_i[i]) for i in H_i])
            n_error_i = np.array([stats.sem(H_i[i]) for i in H_i])
            N = sum(n_i)
            # print(n_i)
            p_i = n_i / N
            p_error_i = n_error_i
            # print(p_i)
            S = N/len(df[df['rho'] == rho]) * np.sum([p * np.log(1 / p) for p in p_i if p != 0])
            N_error = stats.sem(N_array)
            S_error = np.sqrt(np.sum([((np.log(p_i[i]) + 1) ** 2) * (p_error_i[i] ** 2) + ((p_i[i] * np.log(p_i[i])) ** 2) * (N_error ** 2) for i in range(len(p_i))]))
            H_links.append(S)
            H_links_err.append(S_error)

            l = rho**(-1/d)

            if horizon == 'Rindler':
                A.append(1 / (l ** (d - 2)))
            elif horizon == 'Dynamic':
                A.append(n_sphere_surfacearea(n=d-2, r=1) / (l ** (d - 2)))

        else:
            slope, intercept, r_value, p_value, std_err = stats.linregress(A, H_links)
            print(f'{round(slope, 3)} +- {round(std_err, 3)}')

            plt.errorbar(A, H_links, fmt='o', yerr=H_links_err, label='Data', capsize = 5, ecolor='black')
            plt.plot(A, [slope * a + intercept for a in A], label='Fit', zorder=0)
            plt.ylabel(r'S')
            if d == 3:
                plt.xlabel(f'$A/l$')
            elif d == 4:
                plt.xlabel(f'$A/l^2$')
            plt.legend()
            plt.tight_layout()
            plt.show()
            # plt.savefig(f'lambda_{d}d_{horizon}.pgf')
            # plt.clf()
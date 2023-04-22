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
    golden_ratio = (5**.5 - 1) / 2

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

def flatten(l):
    return list(itertools.chain.from_iterable(l))

horizon_types = ['Rindler', 'Dynamic']
d_array = [2, 3, 4]
for horizon in horizon_types:
    print(horizon)
    for d in d_array:
        df = pd.read_csv(f'H_{horizon}{d}d_lambda.csv', names=['rho', 'H', 'b'], header=None)

        rho_array = df['rho'].unique()
        rho_array.sort()

        H_links = []
        H_links_err = []
        A = []
        links_i = []

        for rho in rho_array:
            links = [count_links(a.replace('[', '').replace(']', '').split(', ')) if a != '[]' else 0 for a in df[df['rho'] == rho]['H'].values]
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
            slope, intercept, r_value, p_value, std_err = stats.linregress(A, H_links)
            print(f'{round(slope, 3)} +- {round(std_err, 3)}')

            plt.errorbar(A, H_links, fmt='o', yerr=H_links_err, label='Data', capsize = 5, ecolor='black')
            plt.plot(A, [slope * a + intercept for a in A], label='Fit', zorder=0)
            plt.ylabel(r'$\langle H_\mathrm{link} \rangle$')
            if d == 4:
                plt.xlabel(f'$A/l^2$')
            plt.legend()
            plt.show()
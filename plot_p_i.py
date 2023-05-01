from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
params = {'font.size': 25,
          'font.family': 'Times New Roman',
          'axes.labelsize': 37,
          'legend.fontsize': 27,
          'xtick.labelsize': 30,
          'ytick.labelsize': 30,
          'figure.figsize': [8.5, 6.5],
          'axes.prop_cycle': plt.cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color']
                                        + ['magenta'])
          }
plt.rcParams.update(params)


bar_index = np.arange(10)
width = 0.35
d_array = [2, 3, 4]
moleculetype = 'lambda'
p_i_Rindler_dict, delta_p_i_Rindler_dict, p_i_Dynamic_dict, delta_p_i_Dynamic_dict = dict(
), dict(), dict(), dict()
for d in d_array:
    p_i = dict()
    for horizon in ['Rindler', 'Dynamic']:
        # df = pd.concat([pd.read_csv(f'H_{BHtype}{d}d_{moleculetype}.csv', names=['rho', 'H'], header=None) for d in d_array])
        df = pd.read_csv(f'H_{horizon}{d}d_lambda_clean.csv', names=[
                         'rho', 'H', 'b'], header=None)
        rho_array = df['rho'].unique()
        rho_array.sort()

        for rho in rho_array:
            H = [a.replace('[', '').replace(']', '').split(
                ', ') if a != '[]' else ['0'] for a in df[df['rho'] == rho]['H'].values]

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

                for i, Hele in enumerate(arrayH):
                    N = np.sum(arrayH)
                    if N == 0:
                        p_i[i+1].append(0)
                    else:
                        p_i[i+1].append(Hele/N)

        std_array = [sem(p_i[i]) for i in p_i.keys()]
        mean_array = [np.mean(p_i[i]) for i in p_i.keys()]

        if horizon == 'Rindler':
            p_i_Rindler_dict[d] = mean_array
            # print(p_i_Rindler_dict)
            delta_p_i_Rindler_dict[d] = std_array
        else:
            p_i_Dynamic_dict[d] = mean_array
            delta_p_i_Dynamic_dict[d] = std_array

        if horizon == 'Rindler':
            padding = [0]*10
        else:
            padding = [width]*10

        # if len(mean_array) >= 4:
        #     plt.bar(bar_index[:5] + padding[:5],
        #             mean_array, width, yerr=std_array, label=horizon, capsize=4)
        #     plt.xticks(bar_index[:5] + width[:5] /
        #                2, ['1', '2', '3', '4', '5'])
        # else:
        #     plt.bar(bar_index[:len(mean_array)] + padding[:len(mean_array)],
        #             mean_array, width, yerr=std_array, label=horizon, capsize=4)
        #     plt.xticks(bar_index[:len(mean_array)] +
        #                width[:len(mean_array)] / 2, ['1', '2', '3', '4'])
    plt.ylabel(r'$p_i$')
    plt.xlabel(r'$i$-lambda')
    plt.legend()
    # plt.tight_layout()
    # plt.savefig(
    #     fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\p_i_against_i_{d}d.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.show()


# %%
params = {'font.size': 25,
          'font.family': 'Times New Roman',
          'axes.labelsize': 37,
          'legend.fontsize': 22,
          'xtick.labelsize': 30,
          'ytick.labelsize': 30,
          'figure.figsize': [8.5, 6.5],
          'axes.prop_cycle': plt.cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color']
                                        + ['magenta'])
          }
plt.rcParams.update(params)


def exponential(x, c, A):
    return A*np.exp(-c*x)


def linear(x, m, c):
    return m*x + c


# Exponential Fit

d_array = [2, 3, 4]
# error propagation dictates dy = dx/x for y = ln(x)

for d in d_array:
    print(f'dimension {d} ')
    # Rindler

    p_i_Rindler = list(p_i_Rindler_dict[d])

    # Remove all trailing 0s
    while p_i_Rindler[-1] == 0:
        p_i_Rindler.pop(-1)

    i = np.arange(1, len(p_i_Rindler)+1)

    delta_p_i_Rindler = list(delta_p_i_Rindler_dict[d])

    while delta_p_i_Rindler[-1] == 0:
        delta_p_i_Rindler.pop(-1)

    popt, pcov = curve_fit(exponential, i, p_i_Rindler,
                           sigma=delta_p_i_Rindler, absolute_sigma=True)
    print(f'c for Rindler is {popt[0]} +- {np.sqrt(pcov[0][0])}')
    xLinspace = np.linspace(1, max(i), 100)
    plt.plot(xLinspace, exponential(xLinspace, *popt),
             label='Rindler Exponential Fit', color='green')
    plt.errorbar(i, p_i_Rindler, yerr=delta_p_i_Rindler, fmt='o', capsize=4,
                 color='#1f77b4', linestyle='', ecolor='#1f77b4', label='Rindler Data')

    # Dynamic
    p_i_Dynamic = list(p_i_Dynamic_dict[d])

    # Remove all trailing 0s
    while p_i_Dynamic[-1] == 0:
        p_i_Dynamic.pop(-1)

    i = np.arange(1, len(p_i_Dynamic)+1)

    delta_p_i_Dynamic = list(delta_p_i_Dynamic_dict[d])

    while delta_p_i_Dynamic[-1] == 0:
        delta_p_i_Dynamic.pop(-1)

    popt, pcov = curve_fit(exponential, i, p_i_Dynamic,
                           sigma=delta_p_i_Dynamic, absolute_sigma=True)
    print(f'c for Dynamic is {popt[0]} +- {np.sqrt(pcov[0][0])}')
    xLinspace = np.linspace(1, max(i), 100)
    plt.plot(xLinspace, exponential(xLinspace, *popt),
             label='Dynamic Exponential Fit', color='red')
    plt.errorbar(i, p_i_Dynamic, yerr=delta_p_i_Dynamic, fmt='o', capsize=4,
                 linestyle='', color='magenta', ecolor='magenta', label='Dynamic Data')

    plt.xticks()
    plt.yticks()
    plt.yscale('log')
    plt.ylabel(r'$\ln(p_i)$')
    plt.xlabel(r'$i$-lambda')
    plt.legend(loc=3)
    plt.savefig(
        fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\p_i_against_i_exponentialfit_{d}d.png', dpi=300, transparent=True, bbox_inches='tight')
    plt.show()
    plt.show()

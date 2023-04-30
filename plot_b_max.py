# import matplotlib
# matplotlib.use('pgf')
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
params = {'font.size': 25,
          'font.family': 'Times New Roman',
          'axes.labelsize': 35,
          'legend.fontsize': 25,
          'xtick.labelsize': 28,
          'ytick.labelsize': 28,
          'figure.figsize': [8.5, 6.5],
          'axes.prop_cycle': plt.cycler(color=plt.rcParams['axes.prop_cycle'].by_key()['color']
                                        + ['magenta'])
          }
plt.rcParams.update(params)

horizon_types = ['Rindler', 'Dynamic']
molecule_types = ['lambda', 'v']
d_array = [2, 3, 4]
for horizon in horizon_types:
    for molecule in molecule_types:
        for d in d_array:
            if molecule == 'lambda':
                df = pd.read_csv(f'H_{horizon}{d}d_lambda_clean.csv', names=[
                                 'rho', 'H', 'b'], header=None)
            elif molecule == 'v':
                df = pd.read_csv(f'H_{horizon}{d}d_v_clean.csv', names=[
                                 'rho', 'H', 'subgraphs', 'connected', 'b'], header=None)

            df['rho'] = df['rho'].round(1)
            dfb = df[df['b'] != 0]
            bList = list(dfb['b'].dropna())
            colors = ['red', 'blue', 'green']
            binsCount = [10, 10, 7]
            # plt.hist(bList, range=(0, 3.5), bins=20, density=True,
            #          histtype='stepfilled', label=d, alpha=0.5)
            plt.hist(bList, bins=binsCount[int(d-2)], density=True, histtype='stepfilled',
                     stacked=True, color=colors[int(d-2)], label=f'{d-1}+1 Data', alpha=0.5)

        plt.ylabel('Normalised Frequency')
        plt.xlabel(r'$b_\mathrm{max}$')
        plt.xticks([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5])
        plt.yticks([1, 2, 3, 4, 5])
        plt.legend()
        plt.tight_layout()
        plt.savefig(
            fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\bmax_{horizon}_{molecule}.png', dpi=300, bbox_inches='tight', transparent=True)
        plt.show()
        plt.clf()

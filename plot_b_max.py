import matplotlib
matplotlib.use('pgf')
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})

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

plt.figure(figsize=set_size(483.69687 * 0.45))

import pandas as pd

horizon_types = ['Rindler', 'Dynamic']
molecule_types = ['lambda', 'v']
d_array = [2,3,4]
for horizon in horizon_types:
    for molecule in molecule_types:
        for d in d_array:
            if molecule == 'lambda':
                df = pd.read_csv(f'H_{horizon}{d}d_lambda.csv', names=['rho', 'H', 'b'], header=None)
            elif molecule == 'v':
                df = pd.read_csv(f'H_{horizon}{d}d_v.csv', names=['rho', 'H', 'subgraphs', 'connected','b'], header=None)

            df['rho'] = df['rho'].round(1)
            dfb = df[df['b'] != 0]
            bList= list(dfb['b'].dropna())
            plt.hist(bList, range=(0,3.5), bins=20, density=True, histtype='stepfilled', label=d, alpha=0.5)
        
        plt.ylabel('Normalised Frequency')
        plt.xlabel(r'$b_\mathrm{max}$')
        plt.legend(title=r'$d$')
        # plt.show()
        plt.savefig(f'b_{horizon}_{molecule}.pgf', format='pgf')
        plt.clf()
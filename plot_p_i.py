import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
plt.rc('font', family='Arial', size=20)

bar_index = np.arange(5)
width = 0.35
d_array = [2]
moleculetype = 'lambda'
p_i_Rindler_dict, delta_p_i_Rindler_dict, p_i_Dynamic_dict, delta_p_i_Dynamic_dict = dict(), dict(), dict(), dict()
for d in d_array:
    for BHtype in ['Rindler', 'Dynamic']:
        #df = pd.concat([pd.read_csv(f'H_{BHtype}{d}d_{moleculetype}.csv', names=['rho', 'H'], header=None) for d in d_array])
        df = pd.read_csv(f'Full_Data/H_{BHtype}{d}d_{moleculetype}.csv', names=['rho', 'H', 'b'], header=None) 
        H_dict = {}
        N = df.shape[0]
        maxrho = max(df['rho'])
        maxrhodf = df[df['rho'] == maxrho]
        
        # Initialise p_i dictionary
        p_i_dict = dict()

        for i in range(10): #sets i_max
            p_i_dict[i+1] = list() 
        
        for index, row in maxrhodf.iterrows():
            H = [int(x) if x != '' else 0 for x in row['H'].replace('[', '').replace(']', '').split(', ')]
    
            # Enter realisation p_i into dictionary
            N = np.sum(H)
            try:
                p_i = np.array(H)/ N
            except: 
                p_i = 0
            for i in range(10):
                try:
                    p_i_dict[i+1].append(p_i[i])
                except: 
                    p_i_dict[i+1].append(0)
                    
        mean_array = []
        std_array = []
        for p_i in p_i_dict.values():
            p_i = np.pad(p_i, (0, len(p_i) - len(p_i)), 'constant')
            p_i = np.nan_to_num(p_i)
            mean_array.append(np.mean(p_i))
            std_array.append(sem(p_i))
        
            if BHtype == 'Rindler':
                p_i_Rindler_dict[d] = mean_array 
                #print(p_i_Rindler_dict)
                delta_p_i_Rindler_dict[d] = std_array
            else: 
                p_i_Dynamic_dict[d] = mean_array
                delta_p_i_Dynamic_dict[d] = std_array
        std_array /= np.sum(mean_array)
        mean_array /= np.sum(mean_array)
        if BHtype == 'Rindler':
            padding = 0
        else:
            padding = width
        plt.bar(bar_index + padding, mean_array[:5], width, yerr=std_array[:5], label=BHtype, capsize = 4)
    plt.ylabel(r'$p_i$')
    plt.xlabel(r'$i$-lambda')
    plt.xticks(bar_index + width / 2, ['1', '2', '3', '4', '5'])
    plt.legend()
    #plt.tight_layout()
    plt.savefig(fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\p_i_against_i_{d}d.png', dpi = 300, transparent = True, bbox_inches='tight')
    plt.show()


#%%

from scipy.optimize import curve_fit
def exponential(x,c,A): 
    return A*np.exp(-c*x)

def linear(x, m, c): 
    return m*x + c


# log linear fit vs exponential fit???
d_array = [2,3,4]
#error propagation dictates dy = dx/x for y = ln(x)

for d in d_array:
    print(f'dimension {d} ')
    # Rindler
    p_i_Rindler = p_i_Rindler_dict[d]
    p_i_Rindler = np.nan_to_num(p_i_Rindler)
    p_i_Rindler = [p_i for p_i in p_i_Rindler if p_i != 0]
    i = np.arange(1, len(p_i_Rindler)+1)
    plt.plot(i, np.log(p_i_Rindler), 'x', label = 'Rindler Data', color = 'red')
    
    delta_p_i_Rindler = delta_p_i_Rindler_dict[d]
    delta_p_i_Rindler = [p_i for p_i in p_i_Rindler if p_i != 0]
    popt, pcov = curve_fit(linear, i, np.log(p_i_Rindler), sigma = np.array(delta_p_i_Rindler)/ np.array(p_i_Rindler), absolute_sigma=True)
    print(f'c for Rindler is {popt[0]} +- {np.sqrt(pcov[0][0])}')
    xLinspace = np.linspace(1,max(i), 100)
    plt.plot(xLinspace, linear(xLinspace, *popt), label = 'Rindler Log-Linear Fit', color = 'orange')
    plt.errorbar(i, np.log(p_i_Rindler),yerr = np.array(delta_p_i_Rindler)/ np.array(p_i_Rindler), capsize = 4, linestyle = '', color = 'red')
    
    
    # Dynamic
    p_i_Dynamic = p_i_Dynamic_dict[d]
    p_i_Dynamic = np.nan_to_num(p_i_Dynamic)
    p_i_Dynamic = [p_i for p_i in p_i_Dynamic if p_i != 0] #Remove 0s to avoid np.log(0) = inf
    i = np.arange(1, len(p_i_Dynamic)+1)
    plt.plot(i, np.log(p_i_Dynamic), 'x', label = 'Dynamic Data', color = 'blue')
    
    delta_p_i_Dynamic = delta_p_i_Dynamic_dict[d]
    delta_p_i_Dynamic= [p_i for p_i in p_i_Dynamic if p_i != 0] 
    popt, pcov = curve_fit(linear, i, np.log( p_i_Dynamic), sigma =  np.array(delta_p_i_Dynamic)/ np.array(p_i_Dynamic), absolute_sigma=True)
    print(f'c for Dynamic is {popt[0]} +- {np.sqrt(pcov[0][0])}')
    xLinspace = np.linspace(1,max(i), 100)
    plt.plot(xLinspace, linear(xLinspace, *popt), label = 'Dynamic Log-Linear Fit', color = 'green')
    plt.errorbar(i, np.log(p_i_Dynamic),yerr =  np.array(delta_p_i_Dynamic)/ np.array(p_i_Dynamic), capsize = 4, linestyle = '', color = 'blue')
    
    plt.xticks(fontsize = 28)
    plt.yticks(fontsize = 20)
    plt.ylabel(r'$\ln(p_i)$', fontsize = 24)
    plt.xlabel(r'$i$-lambda', fontsize = 24)
    plt.legend(fontsize = 14, loc = 1)
    plt.savefig(fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\p_i_against_i_exponentialfit_{d}d.png', dpi = 300, transparent = True, bbox_inches='tight')
    plt.show()
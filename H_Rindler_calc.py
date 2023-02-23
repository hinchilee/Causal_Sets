import pandas as pd
import itertools
from causalsetfunctions import find_entropy, linear
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

d_array = [2,3,4]
for d in d_array: 
    df = pd.read_csv(f'H_Rindler{d}d.csv', names=['rho', 'H'], header=None)
    #print(df)
    rho_array = df['rho'].unique()
    rho_array.sort()
    y_entropyList = list() 
    x = list()
    for rho in rho_array:
        
        print('\n sprinkling density',rho, f'in {d} dimensions')
        iterationsNo = df[df['rho'] == rho].shape[0]
        print(f'number of iterations: {iterationsNo}')
        totalHarray = [sum(x) for x in itertools.zip_longest(*[[int(x) if x != '' else 0 for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)]
        print('total Harray \n', totalHarray)
        print('p_i:', totalHarray/ np.sum(totalHarray))
        
        totalLinks = 0 
        for i, value in enumerate(totalHarray):
            totalLinks += (i+1)*value
            
        print(f'Total Links: {totalLinks}')
        
        empiricalavalue = rho**((2-d)/d)*(totalLinks/iterationsNo) 
        
        print(f'Empirical a value {empiricalavalue}')
        
        entropy = find_entropy(totalHarray, iterationsNo)
        print(f'Entropy: {entropy}')
         
        
        y_entropyList.append(entropy)  #S_boltzmann agaisnt A/rho*
        #y_entropyList.append(sum(totalHarray)/iterationsNo) #<N> against A/rho*
        x.append(1/rho**((2-d)/d))
    
    plt.scatter(x, y_entropyList)
    plt.xlabel(r'$\frac{A}{\mathcal{l} ^{ d-2}}$', fontsize = 20 )
    plt.ylabel(r'$s_{Boltz}$', fontsize = 20 )
    popt, pcov = curve_fit(linear, x, y_entropyList)
    xLinspace = np.linspace(min(x), max(x), 100)
    plt.plot(xLinspace, linear(xLinspace, *popt), label = f'{d}d')
    #plt.title(f's_Boltzmann in Rindler in {d}d')
    
    print(f'\n \n \n a_Boltzmann value for {d}d is {popt[0]}')

plt.legend()    
plt.show() 

        
        
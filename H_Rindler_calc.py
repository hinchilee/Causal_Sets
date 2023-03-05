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
        dataArray = np.array([x for x in itertools.zip_longest(*[[int(x) if x != '' else 0 for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)])
        print('total Harray \n', totalHarray)
        print('p_i:', totalHarray/ np.sum(totalHarray))
        
        totalLinks = 0 
        for i, value in enumerate(totalHarray):
            totalLinks += (i+1)*value
            
        print(f'Total Links: {totalLinks}')
        
        empiricalavalue = rho**((2-d)/d)*(totalLinks/iterationsNo)
        
        dataArrayLinks = dataArray
        for row in range(dataArray.shape[0]): 
            dataArrayLinks[row,:] = dataArray[row,:]*(row+1)
        LinksArray = np.sum(dataArrayLinks, axis = 0)
        percaErr = np.std(LinksArray)/ (totalLinks/iterationsNo)
        #due to flucutations in std<H1>
        aerror = percaErr*empiricalavalue/ np.sqrt(iterationsNo)
        print(f'Empirical a value {empiricalavalue} +- {aerror} ')
        
        entropy = find_entropy(totalHarray, iterationsNo)
        #due to fluctiations in <N>, avr no. of molecules per realisation
        MoleculeArray = np.sum(dataArray, axis = 0) 
        percEntropyError = np.std(MoleculeArray)/ ((np.sum(totalHarray)/ iterationsNo))
        entropyerror = percEntropyError * entropy/ np.sqrt(iterationsNo)
        print(f'Entropy: {entropy} +- {entropyerror}')
         
        
        y_entropyList.append(entropy)  #S_boltzmann agaisnt A/rho*
        #y_entropyList.append(sum(totalHarray)/iterationsNo) #<N> against A/rho*
        x.append(1/rho**((2-d)/d))
    if d==4:
        plt.rc('font', family='Arial')
        plt.scatter(np.array(x), np.array(y_entropyList), label = 'Data') 
        plt.errorbar(np.array(x), np.array(y_entropyList), yerr = entropyerror, capsize = 4, linestyle = '')
        plt.xlabel(r'$A/\ell^{d-2}$', fontsize = 25)
        plt.ylabel(r'$s_{Boltz}$', fontsize = 25 )
        popt, pcov = curve_fit(linear, np.array(x), np.array(y_entropyList))
        xLinspace = np.linspace(min(np.array(x)), max(np.array(x)), 100)
        plt.plot(xLinspace, linear(xLinspace, *popt), label = 'Linear Fit', color = 'red')
        #plt.title(f's_Boltzmann in Rindler in {d}d')
    
        print(f'\n \n \n a_Boltzmann value for {d}d is {popt[0]} +- {np.sqrt(pcov[0][0])}')
plt.title('Boltzmannian Entropy for 3+1 Rindler', fontsize = 25, pad = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.legend(fontsize = 15)    
plt.savefig('BoltzEntropyRindler.png', dpi = 300, bbox_inches='tight')
plt.show() 

        
        
import pandas as pd
import itertools
from causalsetfunctions import find_entropy, linear
import numpy as np
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

d_array = [2,3,4]
moleculetype = 'v'
for d in d_array: 
    if moleculetype == 'lambda':
        df = pd.read_csv(f'H_Rindler{d}d_lambda.csv', names=['rho', 'H'], header=None)
        df['rho'] = df['rho'].round(1)  
    elif moleculetype == 'v':
        df = pd.read_csv(f'H_Rindler{d}d_v.csv', names=['rho', 'H', 'subgraphs', 'connected'], header=None)
        df['rho'] = df['rho'].round(1)  
    #print(df)
    rho_array = df['rho'].unique()
    rho_array.sort()
    y_entropyList = list() 
    x = list()
    for rho in rho_array:
        
        print('\n sprinkling density',rho, f'in {d} dimensions')
        iterationsNo = df[df['rho'] == rho].shape[0]
        print(f'number of iterations: {iterationsNo}')
        if moleculetype == 'lambda':
            totalHarray = [sum(x) for x in itertools.zip_longest(*[[int(x) if x != '' else 0 for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)]
            dataArray = np.array([x for x in itertools.zip_longest(*[[int(x) if x != '' else 0 for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)])
        elif moleculetype == 'v':
            totalHarray = np.sum([df[df['rho'] == rho]['H'].values])
            dataArray = [df[df['rho'] == rho]['H'].values]
        print('total Harray \n', totalHarray)
        print('p_i:', totalHarray/ np.sum(totalHarray))
        if moleculetype == 'lambda':
            totalLinks = 0 
            for i, value in enumerate(totalHarray):
                totalLinks += (i+1)*value
                
            print(f'Total Links: {totalLinks}')
    
        elif moleculetype == 'v': 
            totalLinks = totalHarray
        
        empiricalavalue = rho**((2-d)/d)*(totalLinks/iterationsNo)
        
        dataArrayLinks = dataArray
        try:
            for row in range(len(dataArray)): 
                dataArrayLinks[row,:] = dataArray[row,:]*(row+1)
        except:
            pass
        LinksArray = np.sum(dataArrayLinks, axis = 0)
        percaErr = np.std(LinksArray)/ (totalLinks/iterationsNo)
        #due to flucutations in std<H1>
        aerror = percaErr*empiricalavalue/ np.sqrt(iterationsNo)
        print(f'Empirical a value {empiricalavalue} +- {aerror} ')
       
        if moleculetype == 'lambda':
            #theoryauncorrected
            if d == 4:
                theoryacorrected =0.173205  
            elif d== 3: 
                theoryacorrected = 0.2188 
            elif d== 2: 
                theoryacorrected = 1/3
            print(f'theoretical a value for rho {d}d is {theoryacorrected} ')
        
        if moleculetype == 'lambda':
            #lambda molecules
            #entropy = find_entropy(totalHarray, iterationsNo)
            #links 
            entropy = totalLinks/iterationsNo
            
        elif moleculetype == 'v': 
            entropy = totalHarray/ iterationsNo
            
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
plt.title(f'Boltzmannian Entropy for 3+1 Rindler for {moleculetype} molecules', fontsize = 25, pad = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.legend(fontsize = 15)    
plt.savefig(f'BoltzEntropyRindler_{moleculetype}.png', dpi = 300, bbox_inches='tight')
plt.show() 

        
        
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 20:07:26 2023

@author: leehi
"""
import numpy as np 
import pandas as pd
import itertools
from causalsetfunctions import n_sphere_surfacearea, find_entropy, linear
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

T = 1
d_array = [2,3,4]
moleculetype = 'lambda'
for d in d_array: 
    if moleculetype == 'lambda':
        df = pd.read_csv(f'H_Dynamic{d}d_lambda.csv', names=['rho', 'H', 'b'], header=None)
        df['rho'] = df['rho'].round(1)  
    #print(df)
    elif moleculetype == 'v':
        df = pd.read_csv(f'H_Dynamic{d}d_v.csv', names=['rho', 'H',  'subgraphs', 'connected', 'b'], header=None)
        df['rho'] = df['rho'].round(1) 
    rho_array = df['rho'].unique()
    rho_array.sort()
    y_entropyList = list() 
    x = list()
    empiricala = list()
    empiricalaerror = list()
    
    totalLinksList = list() 
    AoverlList = list()
    totalLinksErrorList = list()
    
    for rho in rho_array:
        
        print('\nsprinkling density',rho, f'in {d} dimensions')
        iterationsNo = df[df['rho'] == rho].shape[0]
        print(f'number of iterations: {iterationsNo}')
        print(f'T = {T}')
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
        
        l = rho**(-1/d)
        modified_surface_area = n_sphere_surfacearea(n = d-2, r = T-1.5*l)
            
        empiricalavalue = rho**((2-d)/d)*(totalLinks/iterationsNo)/ n_sphere_surfacearea(n = d - 2, r = T)
        print('surface area:', n_sphere_surfacearea(n = d-2, r= T))

        totalLinksList.append(totalLinks/iterationsNo)
        AoverlList.append((n_sphere_surfacearea(n = d-2, r= T)/ rho**((2-d)/d)))
        
        if d== 4:
            empiricala.append(empiricalavalue)
        
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
        if d == 4:
            empiricalaerror.append(aerror)
        print(f'Empirical a value {empiricalavalue} +- {aerror} ')
        
        totalLinksErrorList.append(np.std(LinksArray)/(totalLinks/iterationsNo)) 
       

        if moleculetype == 'lambda':
            #theoryauncorrected
            if d == 4:
                theoryacorrected =0.173205 + (-0.0209261)*(d-2)*(rho**(-1/d))/(T*np.sqrt(2)) 
            elif d== 3: 
                theoryacorrected = 0.2188 
            elif d== 2: 
                theoryacorrected = 2/3
            print(f'theoretical a value for rho {d}d after finite rho correction is {theoryacorrected} ')
    
        if moleculetype == 'lambda':
            #lambda molecules
            #entropy = find_entropy(totalHarray, iterationsNo)
            #links 
            entropy = totalLinks/iterationsNo
        elif moleculetype == 'v': 
            entropy = totalHarray/iterationsNo
        #due to fluctiations in <N>, avr no. of molecules per realisation
        MoleculeArray = np.sum(dataArray, axis = 0) 
        percEntropyError = np.std(MoleculeArray)/ ((np.sum(totalHarray)/ iterationsNo))
            
        entropyerror = percEntropyError * entropy/ np.sqrt(iterationsNo)
        print(f'Entropy: {entropy} +- {entropyerror}')
        
        y_entropyList.append(entropy)
        #y_entropyList.append(sum(totalHarray)/iterationsNo) #<N> against A/rho*
        x.append(n_sphere_surfacearea(n = d - 2, r = T)/rho**((2-d)/d))
    
    if moleculetype == 'lambda':
        # Plots the link molecules of <H_links> against A/rho**(2-d/d)    
        plt.scatter(AoverlList, totalLinksList, label = 'Data')
        plt.errorbar(AoverlList, totalLinksList, yerr = totalLinksErrorList, capsize = 4, linestyle = '')
        popt, pcov = curve_fit(linear, AoverlList, totalLinksList)
        xLinspace = np.linspace(min(AoverlList), max(AoverlList), 100)
        plt.plot(xLinspace, linear(xLinspace, *popt), label = 'Linear Fit', color = 'red') 
        plt.xlabel(r'$A/{\l^{d-2}}$', fontsize = 25)
        plt.ylabel(r'$\langle H_{links} \rangle$', fontsize = 25 )
        print(f'\n \n \n a_Boltzmann value for {d}d is {popt[0]} +- {np.sqrt(pcov[0][0])}')
    
        #plt.title(f'Link Counting for {d-1}+1 Dynamic', fontsize = 25, pad = 20)
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        plt.legend(fontsize = 15)    
        plt.savefig(fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\LinksEntropyDynamic_{moleculetype}_{d}d.png', dpi = 300, bbox_inches='tight', transparent = True)
        plt.show() 
        
        
        # plot v-molecules or lambda molecules
    plt.rc('font', family='Arial')
    #x = x[:-1]
    #y_entropyList = y_entropyList[:-1]
    plt.scatter(np.array(x), np.array(y_entropyList), label = 'Data') 
    plt.errorbar(np.array(x), np.array(y_entropyList), yerr = entropyerror, capsize = 4, linestyle = '')
    plt.xlabel(r'$A/\ell^{d-2}$', fontsize = 25)
    plt.ylabel(r'$s_{Boltz}$', fontsize = 25 )
    #plt.title(f'Boltzmannian Entropy for {d-1}+1 Dynamic', fontsize = 25, pad = 20)
    if moleculetype == 'v':
        plt.ylabel(r'$\langle H_V \rangle$', fontsize = 25 )
    elif moleculetype == 'lambda': 
        plt.ylabel(r'$S_{Boltz}$', fontsize = 25 )
    popt, pcov = curve_fit(linear, np.array(x), np.array(y_entropyList))
    xLinspace = np.linspace(min(np.array(x)), max(np.array(x)), 100)
    plt.plot(xLinspace, linear(xLinspace, *popt), label = 'Linear Fit', color = 'red')
    #plt.title(f's_Boltzmann in Rindler in {d}d')
    
    print(f'\n \n \n a_Boltzmann value for {d}d is {popt[0]} +- {np.sqrt(pcov[0][0])}')
    
    
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)
    plt.legend(fontsize = 15)  
    plt.savefig(fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\BoltzEntropyDynamic_{moleculetype}_d{d}.png', dpi = 300, bbox_inches='tight', transparent = True)  
    plt.show() 
        
        
 #%%       
for d in d_array: 
    if moleculetype == 'lambda':
        df = pd.read_csv(f'H_Dynamic{d}d_lambda.csv', names=['rho', 'H', 'b'], header=None)
        df['rho'] = df['rho'].round(1)  
    #print(df)
    elif moleculetype == 'v':
        df = pd.read_csv(f'H_Dynamic{d}d_v.csv', names=['rho', 'H',  'subgraphs', 'connected', 'b'], header=None)
        df['rho'] = df['rho'].round(1) 
    #   Analyse b 
    dfb = df[df['b'] != 0] #dropped 0 (optional)
    bList= list(dfb['b'].dropna()) #dropped Nans
    colors = ['red', 'blue', 'green']
    plt.hist(bList, bins = 10, density = True, histtype = 'stepfilled', stacked = True, color = colors[int(d-2)], label = f'{d-1}+1 Data', alpha = 0.5)

plt.ylabel('Normalised Frequency')
plt.xlabel(r'$b$')
plt.legend()
plt.savefig(fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\bepislonDistribution_Dynamic_{moleculetype}.png', dpi = 300, bbox_inches='tight', transparent = True)
plt.show()
        

#%%

plt.rc('font', family='Arial')
plt.scatter(rho_array, empiricala, label = 'Data') 
plt.errorbar(rho_array, empiricala, yerr = empiricalaerror, capsize = 4, linestyle = '')
plt.xlabel(r'$\rho$', fontsize = 25)
plt.ylabel(r'$a_{empiricial}$', fontsize = 25 )
#popt, pcov = curve_fit(linear, np.array(x), np.array(y_entropyList))
#xLinspace = np.linspace(min(np.array(x)), max(np.array(x)), 100)
#plt.plot(xLinspace, linear(xLinspace, *popt), label = 'Linear Fit', color = 'red')
#plt.title(f's_Boltzmann in Rindler in {d}d')

#print(f'\n \n \n a_Boltzmann value for {d}d is {popt[0]} +- {np.sqrt(pcov[0][0])}')

plt.title(r'Empirical $a^{(4)}$ for 3+1 Dynamic', fontsize = 25, pad = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.legend(fontsize = 25, loc = 4)  
plt.savefig(fr'C:\Users\leehi\OneDrive\Documents\Imperial_tings\Fourth_Year\MSci Project\Thesis\Plots\Fittedvalue_a4_{d}D_Dynamic.png', dpi = 300, bbox_inches='tight', transparent = True)  
plt.show() 


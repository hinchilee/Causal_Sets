import pandas as pd
import itertools

d_array = [3]
for d in d_array: 
    df = pd.read_csv(f'H_Rindler{d}d.csv', names=['rho', 'H'], header=None)
    #print(df)
    rho_array = df['rho'].unique()
    for rho in rho_array:
        
        print('\n sprinkling density',rho, f'in {d} dimensions')
        iterationsNo = df[df['rho'] == rho].shape[0]
        print(f'number of iterations: {iterationsNo}')
        totalHarray = [sum(x) for x in itertools.zip_longest(*[[int(x) for x in a.replace('[', '').replace(']', '').split(', ')] for a in df[df['rho'] == rho]['H'].values], fillvalue=0)]
        print('total Harray \n', totalHarray)
        
        totalLinks = 0 
        for i, value in enumerate(totalHarray):
            totalLinks += (i+1)*value
            
        print(f'Total Links: {totalLinks}')
        
        empiricalavalue = rho**((2-d)/d)*(totalLinks/iterationsNo) 
        
        print(f'Empirical a value {empiricalavalue}')
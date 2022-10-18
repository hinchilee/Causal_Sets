
#Causal Set Methods

from scipy.optimize import fsolve
from scipy.special import gamma 

def MM_dimension(d, r): 
    #Solves d_MM(r) = 0
    
    if r < 0 or r > 1: 
        raise ValueError('Make sure ordering fraction, r is between 0 and 1!')
        
    return 1.5*gamma(d/2)*gamma(d)/ gamma(3*d/2) - r

root = fsolve(MM_dimension, x0 = 4, args = (0.05))
print(root)

def findOrderingFraction(Interval, CausalMatrix):
    
    #get element labels 
    #check pairwise (itertools pair??) 
    #for each pair, 
        #add value of element in CausalMatrix to running count 
    #final count = 2R
    
    n = len(Interval)
    r = 2*R/ (n*(n-1))
    return r 

def to_the_past(ele1, ele2, CausalMatrix): 
    '''Checks whether ele1 < ele2'''
    
    l1 = ele1.label 
    l2 = ele2.label 
    
    if CausalMatrix[l1][l2] == 1:
        return True
    else: 
        return False
    
def to_the_future(ele1, ele2, CausalMatrix): 
    '''Checks whether ele1 > ele2'''
    
    l1 = ele1.label 
    l2 = ele2.label 
    
    if CausalMatrix[l2][l1] == 1:
        return True
    else: 
        return False
    
def spacelike_to(ele1, ele2, CausalMatrix): 
    '''Checks whether ele1 is spacelike to ele2'''
    
    l1 = ele1.label 
    l2 = ele2.label 
    
    if CausalMatrix[l1][l2] == 0 and CausalMatrix[l2][l1] == 0:
        return True
    else: 
        return False
    
    
def check_relation(ele1, ele2, CausalMatrix): 
    '''gives relation R in ele1 R ele2, where R = {precedes, succeeds, spacelike to}'''
    
    l1 = ele1.label 
    l2 = ele2.label 
    if CausalMatrix[l1][l2] == 1:
        state = 'precedes'
    elif CausalMatrix[l2][l1] == 1: 
        state = 'succedes'
    else: 
        state = 'is spacelike to'
    print(f'Causal Event {l1} {state} Causal Event {l2}.')
    
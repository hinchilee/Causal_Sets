
import numpy as np 
from math import gamma

#Causal Set Methods

def spacetime_interval(stcoord1, stcoord2, BHtype , wrapAroundLength = 2): 

    #Returns the spacetime interval given two spacetime coordinates
    diff = (stcoord1 - stcoord2)**2
    spacetime_interval = -diff[0] #time interval 
    
    if BHtype == 'Rindler': 
        #Hardcoded to wrap around y, z
        for spacedimension in range(1, len(diff)): 
            if spacedimension == 2 or spacedimension == 3:
                spacetime_interval += min(diff[spacedimension], (wrapAroundLength - np.sqrt(diff[spacedimension]))**2)
            else:
                spacetime_interval += diff[spacedimension]
    
    elif BHtype == 'Dynamic':
        #Hardcoded to wrap around y 
        for spacedimension in range(1, len(diff)): 
            if spacedimension == 2:
                spacetime_interval += min(diff[spacedimension], (wrapAroundLength - np.sqrt(diff[spacedimension]))**2)
            else:
                spacetime_interval += diff[spacedimension]
    
    elif BHtype == 'Empty': 
        for spacedimension in range(1, len(diff)): 
            spacetime_interval += diff[spacedimension]
    
    return spacetime_interval

def find_entropy(H_array):
    N = np.sum(H_array)
    p_array = H_array / N
    
    return N * np.sum([p_i * np.log(1 / p_i) for p_i in p_array if p_i != 0])

def inside_horizon(x):
    s = -(x[0] ** 2)
    for x_i in x[1:]:
        s += x_i ** 2

    if s < 0:
        return True
    elif s > 0:
        return False
    
def convert_Harray_1molecule(H_array): 
	ans = 0 
	for i, count in enumerate(H_array): 
		ans += i*count
	return ans

def find_entropy_1molecule(oneMoleculeCount): 
    return np.log(oneMoleculeCount)

if __name__ == "__main__":
    H_arr = np.array([2,1])
    print(find_entropy(H_arr))
    # c = np.array([[0,0,0,1,1],
    #               [0,0,0,1,1],
    #               [0,0,0,1,1],
    #               [0,0,0,0,1],
    #               [0,0,0,0,0]])
    # LinkMatrix = c - fC2(c)
    # print(LinkMatrix)
    #print(spacetime_interval(np.array([1,1.1]), np.array([0,0])))
    


import numpy as np 
from math import gamma

#Causal Set Methods

def spacetime_interval(stcoord1, stcoord2, periodicBC = False, wrapAroundLength = 2): 
    wraparound_dimension = 2

    #Returns the spacetime interval given two spacetime coordinates
    diff = (stcoord1 - stcoord2)**2
    spacetime_interval = -diff[0] #time interval 
    for spacedimension in range(1, len(diff)): 
        if periodicBC == True: 
            if spacedimension == wraparound_dimension:
                spacetime_interval += min(diff[spacedimension], (wrapAroundLength - np.sqrt(diff[spacedimension]))**2)
            else:
                spacetime_interval += diff[spacedimension]
        else: 
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
    


import numpy as np 

#Causal Set Methods

def fslow(ele): 
    # f(x) = 0 if x = 0 ; 1 otherwise
    if ele == 0: 
        return 0
    else: 
        return 1

f = np.vectorize(fslow)
    
def findC2(CausalMatrix):
    #Intermediate step for finding Link Matrix
    #Matrix that counts the number of 3-chains between elements i j
    #Optimised to consider only top half of matrix
    
    C2: np.ndarray = np.zeros(CausalMatrix.shape)
    for i in range(len(CausalMatrix)):
        for j in range(i,len(CausalMatrix)):
            count = 0 
            for k in range(len(CausalMatrix)): 
                count += CausalMatrix[i, k]*CausalMatrix[k, j]
            C2[i,j] = count  
    return C2 
    
def spacetime_interval(stcoord1, stcoord2): 
    #Returns the spacetime interval given two spacetime coordinates
    
    diff = (stcoord1 - stcoord2)**2
    spacetime_interval = -diff[0] #time interval 
    for spacedimension in range(1, len(diff)): 
        spacetime_interval += diff[spacedimension]
    return spacetime_interval

if __name__ == "__main__":
    print(spacetime_interval(np.array([1,1.1]), np.array([0,0])))
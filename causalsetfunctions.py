
import numpy as np 

#Causal Set Methods
    
def fC2(CausalMatrix):
    # f(x) = 0 if x = 0 ; 1 otherwise
    #Intermediate step for finding Link Matrix
    #Matrix that counts the number of 3-chains between elements i j
    #Optimised to consider only top half of matrix
    
    C2: np.ndarray = np.zeros(CausalMatrix.shape)
    for i in range(len(CausalMatrix)):
        for j in range(i,len(CausalMatrix)):
            count = np.sum(CausalMatrix[i, :] * CausalMatrix[:, j])
            C2[i,j] = 1 if count >= 1 else 0  
    return C2 

def spacetime_interval(stcoord1, stcoord2): 
    wraparound_dimension = 2

    #Returns the spacetime interval given two spacetime coordinates
    diff = (stcoord1 - stcoord2)**2
    spacetime_interval = -diff[0] #time interval 
    for spacedimension in range(1, len(diff)): 
        if spacedimension == wraparound_dimension:
            spacetime_interval += min(diff[spacedimension], (2 - np.sqrt(diff[spacedimension]))**2)
        else:
            spacetime_interval += diff[spacedimension]
    return spacetime_interval

def generate_CausalMatrix(ElementList): 
    # Induce causal relations by transitivity
    
    A = np.zeros((len(ElementList), len(ElementList)), dtype = int)
    for j in range(len(ElementList)): 
        for i in reversed(range(j)): 
            if A[i,j] == 0:
                if spacetime_interval(ElementList[j], ElementList[i]) < 0: 
                    A[i,j] = 1 
                    #Then inherit i's past 
                    A[:,j] = np.bitwise_or(A[:,j], A[:,i])
                else: 
                    pass 
            else: 
                pass
    return A 

if __name__ == "__main__":
    c = np.array([[0,0,0,1,1],
                  [0,0,0,1,1],
                  [0,0,0,1,1],
                  [0,0,0,0,1],
                  [0,0,0,0,0]])
    LinkMatrix = c - fC2(c)
    print(LinkMatrix)
    #print(spacetime_interval(np.array([1,1.1]), np.array([0,0])))
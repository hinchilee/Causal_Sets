
import numpy as np 

C = np.array([[0,0,0,1,1],
              [0,0,0,1,1],
              [0,0,0,0,1],
              [0,0,0,0,1],
              [0,0,0,0,0]])

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
    
    C2: np.ndarray = np.zeros(CausalMatrix.shape)
    for i in range(len(CausalMatrix)):
        for j in range(len(CausalMatrix)):
            count = 0 
            for k in range(len(CausalMatrix)): 
                count += CausalMatrix[i, k]*CausalMatrix[k, j]
            C2[i,j] = count  
    return C2 

def findLinkMatrix(CausalMatrix):
    #L_ij = 1 if e_i <* e_j; 0 otherwise, where <* is a link
    #L = C - f(C2)
    
    LinkMatrix: np.array = CausalMatrix - f(findC2(CausalMatrix))
    return LinkMatrix

L = findLinkMatrix(C)
print(L)


import numpy as np

def Sprinkling_Minkowski(d = 2, n = 10, bounds = np.array([[0,1],[0,1]])): 
    #Generate a list of n points of (t, x1, x2, ..., x_d) randomly in the bounds
    #bounds is a nparray of (low, high) bounds, first ele is time, the remaining is space
    
    sprinkledCoords:np.array = np.random.rand(n,d)
    
    #Scaling the distribution to fit the time bounds
    sprinkledCoords[:,0] *= (bounds[0][1] - bounds[0][0])
    sprinkledCoords[:,0] += (bounds[0][0])
    
    #Scaling the distribution to fit the space bounds
    for xn in range(1,d):
        sprinkledCoords[:,xn] *= (bounds[xn][1] - bounds[xn][0])
        sprinkledCoords[:,xn] += (bounds[xn][0])
    
    return sprinkledCoords

if __name__ == "__main__":
    s = Sprinkling_Minkowski()
    print(s)
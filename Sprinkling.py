
import numpy as np

def Sprinkling_Uniform(dimension = 2, number_of_points = 10, bounds = np.array([[0,1],[0,1]])): 
    #Generate a list of n points of (t, x1, x2, ..., x_d) randomly in the bounds
    #bounds is a nparray of (low, high) bounds, first ele is time, the remaining is space
    
    if len(bounds) != dimension: 
        raise ValueError('The dimension of bounds does not match the bounds of the sprinkling dimension!')
    
    sprinkledCoords:np.array = np.random.rand(number_of_points,dimension)
    
    #Scaling the distribution to fit the time bounds
    sprinkledCoords[:,0] *= (bounds[0][1] - bounds[0][0])
    sprinkledCoords[:,0] += (bounds[0][0])
    
    #Scaling the distribution to fit the space bounds
    for xn in range(1,dimension):
        sprinkledCoords[:,xn] *= (bounds[xn][1] - bounds[xn][0])
        sprinkledCoords[:,xn] += (bounds[xn][0])
    
    return sprinkledCoords

if __name__ == "__main__":
    s = Sprinkling_Uniforms()
    print(s)

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

def Sprinkling_Bicone(dimension = 2, number_of_points = 10): 
    
    Time = list() 
    Space = list() 
    T = 1 
    R_0 = 1
    for i in range(number_of_points): 
        x = np.random.uniform(0,1)
        sign = np.random.uniform (0,1)
        t = T*(1-x**(1/dimension))
        if sign > 0.5: 
            Time.append(t)
        else: 
            Time.append(-t)
        R_t = R_0* x ** (1/dimension)
        v_unnormalised = np.random.normal(0, 1, size = dimension - 1)
        
        def normalise(v): 
            norm = np.linalg.norm(v)
            if norm == 0: 
                return v 
            return v/norm
        
        v = normalise(v_unnormalised)
        k = np.random.uniform(0,1)
        v_f = v*R_t*k**(1/(dimension - 1))
        Space.append(v_f)
        
    Time = np.array(Time)
    Space = np.array(Space)
    Coords = np.concatenate((Time.reshape(number_of_points, 1), Space), axis = 1)

    #Manually add top and bottom points  
    top_point = np.concatenate(([1],np.zeros(shape = (dimension-1, ), dtype = int)))
    bottom_point = np.concatenate(([-1],np.zeros(shape = (dimension-1, ), dtype = int)))
    Coords = np.concatenate((Coords, top_point.reshape(1,dimension)), axis = 0)
    Coords = np.concatenate((Coords, bottom_point.reshape(1,dimension)), axis = 0)
    
    return Coords

def Sprinkling_Tube(dimension = 2, number_of_points = 10, bounds = [ 0.3, 1, 0, 1] ): 
    
    R_min, R_max, T_min, T_max = bounds
    Time = list() 
    Space = list() 
    x_min = (R_min/ R_max)**(dimension-1) 
    for i in range(number_of_points): 
        x = np.random.uniform(0,1)
        xprime = (1-x_min)*x + x_min
        t2 = np.random.uniform(0,1)
        t = (T_max-T_min)*t2 + T_min
        R_t = R_max* xprime ** (1/(dimension-1))
        v_unnormalised = np.random.normal(0, 1, size = dimension - 1)
        
        def normalise(v): 
            norm = np.linalg.norm(v)
            if norm == 0: 
                return v 
            return v/norm
        
        v = normalise(v_unnormalised)
        k = np.random.uniform(0,1)
        v_f = v*R_t#*k**(1/(dimension - 1))
        #print(k**(1/(dimension - 1)))
        Space.append(v_f)
        Time.append(t)
        
    Time = np.array(Time)
    Space = np.array(Space)
    Coords = np.concatenate((Time.reshape(number_of_points, 1), Space), axis = 1)    
    #print(Coords)
    return Coords
    


if __name__ == "__main__":
    np.random.seed(11)
    a= Sprinkling_Tube(dimension = 3, 
                          bounds = np.array([[-10,10],
                                             [-10,10], 
                                             [-10,10]]))
    print(a)

import numpy as np 
from math import gamma
from scipy.special import hyp2f1
import networkx as nx

#Causal Set Methods

def spacetime_interval(stcoord1, stcoord2, BHtype , wrapAroundLength = 2): 

    #Returns the spacetime interval given two spacetime coordinates
    diff = (stcoord2 - stcoord1)*(stcoord2 - stcoord1)
    spacetime_interval = -diff[0] #time interval 
    
    if BHtype == 'Rindler': 
        #Hardcoded to wrap around y, z
        for spacedimension in range(1, len(diff)): 
            if spacedimension == 2 or spacedimension == 3:
                a =  (wrapAroundLength - stcoord2[spacedimension] + stcoord1[spacedimension])
                spacetime_interval += min(diff[spacedimension], a*a)
            else:
                spacetime_interval += diff[spacedimension]
    
    elif BHtype == 'Dynamic':
        #Hardcoded to no wrap around
        for spacedimension in range(1, len(diff)): 
            spacetime_interval += diff[spacedimension]
    
    elif BHtype == 'Empty': 
        for spacedimension in range(1, len(diff)): 
            spacetime_interval += diff[spacedimension]
    
    return spacetime_interval

def find_entropy(H_array, realisations):
    N = np.sum(H_array)
    p_array = H_array / N
    
    return (N/realisations) * np.sum([p_i * np.log(1 / p_i) for p_i in p_array if p_i != 0])

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
    return 2*np.pi*oneMoleculeCount

def theoretical_a4_flat(n):
    # Eqn 3.12 in Barton et al 2019 
    return 2**(8*n+1)/ (3**(3*n+0.5)*np.sqrt(np.pi)*np.math.factorial(n))*gamma(n+0.5)*(B(3/8, 3*n+1, n+1) + B(5/8, n+1, 3*n+1) - B(1/4, n+1, 3*n+1))

def B(z, a, b):
    # Eqn 3.13 in Barton et al 2019
    return (z**a)*hyp2f1(a, 1-b, a+1, z)/a

def frustum_curved_surfaceaarea(h1, h2):
    # Returns curved surface area of frustum, used to calculate A_BH for dynamic black holes 
    # convention h1 < h2
    return np.pi*(h2*(h2+np.sqrt(2*h2**2))- h1*(h1+np.sqrt(2*h1**2)))
    
def two_torus_surfacearea(R, r):
    return 4*np.pi**2*r*R

def n_ball_volume(n, r):
    if n == 0: 
        return 1 
    elif n == 1: 
        return 2*r
    else: 
        return (2*r*r*np.pi/n)*n_ball_volume(n-2, r)
    
def n_sphere_surfacearea(n, r): 
    try:
        ans = 2*(np.pi**((n+1)/2))/(gamma((n+1)/2))
    except: 
        ans = 1
    if n == 0 and r == 1: 
        ans = 2
        
    return ans
    
def linear(x, m, c): 
    return m*x+c

def generate_adjacency_matrix(VElementsLabelsList, CausalMatrix): 
    #AdjMatrix symmetric for undirected graph, 0 on the diagonals
    
    AdjMatrix = CausalMatrix[VElementsLabelsList,:][:,VElementsLabelsList]
    AdjMatrix += AdjMatrix.T
    return AdjMatrix

def count_subgraphs(AdjMatrix):
    G = nx.from_numpy_matrix(AdjMatrix)
    return nx.number_connected_components(G)

def check_connected_graph(AdjMatrix):
    G = nx.from_numpy_matrix(AdjMatrix)
    return nx.is_connected(G)

# def compute_spacetimecuts_uniform_Rindler(d, rho, N_max = 10000, b = 4, bhEdge = 0.5): 
#     #outputs bounds array of dimension d
    
#     l = rho**(-1/d)
#     print('l:', l)
#     epsilon = b*l 
#     stV = (2*epsilon)*epsilon   #wrap around dimensions length = 2epsilon
#     print(stV)

#     boundsArray = np.array([[-bhEdge, bhEdge] for i in range(d)])
#     boundsArray[0][0] = bhEdge - epsilon 
#     if epsilon > 1: 
#         boundsArray[1][1] = bhEdge + bhEdge*2
#     else: 
#         boundsArray[1][1] = bhEdge + epsilon
#         boundsArray[1][0] = bhEdge - epsilon
        
#     stV = np.prod(boundsArray[:,1] - boundsArray[:,0])
#     print(stV)
#     adjusted_rho = N_max/ stV
#     adjusted_l = adjusted_rho**(-1/d)
    
#     # Ensures that b is at least larger than initial value
#     if adjusted_rho < rho: 
#         raise ValueError(f'Adjusted_rho = {adjusted_rho} < rho = {rho}. Please choose a SMALLER rho.')
    
#     return boundsArray, adjusted_rho, l, adjusted_l

def compute_spacetimecuts_uniform_Rindler(d, rho2, N_max = 10000, b = 4, bhEdge = 0.5): 
    #outputs bounds array of dimension d
    
    
    if d == 2: 
        l = rho2**(-1/d)
        b = (N_max/2)**0.5
    elif d == 3: 
        l = 2*b**2/N_max
    elif d == 4: 
        l = (2*b**2/N_max)**0.5
        
    #print('l:', l)
    epsilon = b*l 

    boundsArray = np.array([[-bhEdge, bhEdge] for i in range(d)])
    boundsArray[0][0] = bhEdge - epsilon 
    boundsArray[1][1] = bhEdge + epsilon
    boundsArray[1][0] = bhEdge - epsilon
    if epsilon > 1: 
        boundsArray[1][1] = bhEdge + bhEdge*2
        boundsArray[0][0] = -bhEdge 
        
    stV = np.prod(boundsArray[:,1] - boundsArray[:,0])
    print(stV)
    adjusted_rho = N_max/stV
    
    return boundsArray, adjusted_rho, l

def compute_spacetimecuts_tube(d, rho2 = 10000, N_max = 10000, b = 3): 
    #outputs bounds array of R_min, R_max, T_min, T_max, #T_max = 1 by default
    
    #assumed T_max = 1
    if d == 2: 
        l = rho2**(-1/d)
        b = (N_max/4)**0.5
    elif d == 3: 
        l = 4*np.pi*b**2/N_max
    elif d == 4: 
        a = 4*np.pi/3
        #care, for too small N_max, may yield complex solution for l 
        l = (6*b**2/(N_max/a - 2*(b**4)))**0.5
        while isinstance(l, complex):
            b -= 0.1
            l = (6*b**2/(N_max/a - 2*(b**4)))**0.5

    T_max = 1 
    T_min = T_max - b*l
    if T_min < 0: 
        T_min = 0
    R_min = T_max - b*l
    if R_min < 0: 
        R_min = 0
    R_max = T_max + b*l
    if R_max > T_max*2: 
        R_max = T_max*2
    
    ndimension = d - 1 
    stV = (T_max-T_min)*(n_ball_volume(ndimension, R_max) - n_ball_volume(ndimension, R_min))
    rho = N_max/stV
    #adjustedl = rho**(-1/d)
    #adjustedb = (R_max - T_max)/ adjustedl
    #print(f'adjustedl:{adjustedl}, adjustedb: {adjustedb}')
    
    return [R_min, R_max, T_min, T_max], rho



if __name__ == "__main__":
    # H_arr = np.array([2,1])
    # print(find_entropy(H_arr))
    # print(spacetime_interval(np.array([1,1.1]), np.array([0,0])))
    #print(theoretical_a4_flat(1))  # Does indeed yield Barton et al 2019 (3.16), should be correct
    #print(frustum_curved_surfaceaarea(3, 5))
    #print(n_sphere_surfacearea(0, 1))
    
    #Dynamic
    dimension = 3
    boundsArray, rho= compute_spacetimecuts_tube(d = dimension, rho2 = 3e19, N_max = 9000, b= 3.7)
    R_min, R_max, T_min, T_max = boundsArray
    ndimension = dimension -1
    stV = (T_max-T_min)*(n_ball_volume(ndimension, R_max) - n_ball_volume(ndimension, R_min))
    points = stV * rho
    print(stV, rho, points, boundsArray)
    
    # #Rindler
    # dimensionList = [2]
    # # Limiting case is in 4d
    # for dimension in dimensionList:
    #     print(f'dimension: {dimension}')
    #     boundsArray, adjusted_rho, l = compute_spacetimecuts_uniform_Rindler(d = dimension, rho2 = 1e20, N_max = 1000, b = 4)
    #     print(boundsArray, adjusted_rho)
    #     stV = np.prod(boundsArray[:,1] - boundsArray[:,0])
    #     print('points:', stV* adjusted_rho)
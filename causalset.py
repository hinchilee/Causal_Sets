# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 12:09:46 2022

@author: leehi
"""

import cProfile, pstats, io
from pstats import SortKey
import numpy as np 
import time 
import matplotlib.pyplot as plt
from scipy.stats import poisson 
from scipy.optimize import fsolve
from scipy.special import gamma 
from causalsetfunctions import fC2, spacetime_interval, inside_horizon
from causalEvent import CausalEvent
from Sprinkling import Sprinkling_Uniform, Sprinkling_Bicone

class CausalSet(object): 
    
    def __init__(self, **kwargs): 
        
        self.dimension = kwargs.get('dimension', 2)
        self.ElementList: list (CausalEvent) = list() 
        self.periodicBC = kwargs.get('periodicBC', False)
        self.dynamic = kwargs.get('DynamicBH')
    
        #Sprinkling and sorting by time coordinate 
        
        # FOR RINDLER HORIZON  
        if not self.dynamic:
            # sprinkledcoords = Sprinkling_Bicone(dimension = self.dimension, number_of_points = kwargs.get('number_of_points', 5))
            # Normalised Sprinkling Volume to 1!!! Important to find out <N> then Poisson, not get N from Poisson then scale to area 
            
            self.SpacetimeVolume = 1
            sprinkledcoords = Sprinkling_Uniform(dimension = kwargs.get('dimension', 2),number_of_points = poisson.rvs(kwargs.get('sprinkling_density')), bounds = kwargs.get('bounds', np.array([[-0.5,0.5] for i in range(kwargs.get('dimension', 2))])))
            self.wrapAroundLength = 1

        # FOR DYNAMIC BLACKHOLE HORIZON
        else:
            # Sprinkle uniformly into a slab to calculate entropy for Sigma_t = T 
            # into bounds of t: [T-1, T+1]; 
            # x, y, z = [-T-2, T+2]; wrap around y
            
            self.T = kwargs.get('T')
            spaceBounds = np.array([-self.T-2, self.T+2])
            timeBounds = np.array([self.T-1, self.T+1])
            self.SpacetimeVolume = ((2*(self.T+2))**(self.dimension - 1))*2
            AveragePoints = self.SpacetimeVolume*(kwargs.get('sprinkling_density'))
            boundsArray = np.concatenate((timeBounds, np.tile(spaceBounds, self.dimension - 1)),axis = 0).reshape(self.dimension, 2)
            sprinkledcoords = Sprinkling_Uniform(dimension = self.dimension, 
                                                number_of_points = poisson.rvs(AveragePoints),
                                                bounds = boundsArray) 
            self.wrapAroundLength = 2*(self.T+2)
                        
        #sort by time
        sprinkledcoords = sprinkledcoords[sprinkledcoords[:, 0].argsort()]
        
        for i, coords in enumerate(sprinkledcoords): 
            self.ElementList.append(CausalEvent(label = i, coordinates = coords))
    
        self.LinkMatrix = None
        
        self.CausalMatrix = self.generate_CausalMatrix()
        
    def generate_CausalMatrix(self): 
        # Induce causal relations by transitivity
        
        A = np.zeros((len(self.ElementList), len(self.ElementList)), dtype = int)
        for j in range(len(self.ElementList)): 
            for i in reversed(range(j)): 
                if A[i,j] == 0:
                    if spacetime_interval(self.ElementList[j].coordinates, self.ElementList[i].coordinates, self.periodicBC, self.wrapAroundLength) < 0: 
                        A[i,j] = 1 
                        #Then inherit i's past 
                        A[:,j] = np.bitwise_or(A[:,j], A[:,i])
                    else: 
                        pass 
                else: 
                    pass
        return A 
    
    
    #VISUALIATION 
    
    def visualisation(self):
        coordinates = np.array([x.coordinates for x in self.ElementList])

        if self.dimension == 2:
            if self.LinkMatrix == None:
                self.find_linkmatrix()

            plt.scatter(coordinates[:, 1], coordinates[:, 0], s = 100)
            for i in range(len(self.LinkMatrix)):
                for j in range(len(self.LinkMatrix[i])):
                    if self.LinkMatrix[i][j] == 1:
                        plt.plot([coordinates[i][1], coordinates[j][1]], [coordinates[i][0], coordinates[j][0]], color = 'green')
            plt.xlabel('Space')
            plt.ylabel('Time')
            plt.axis('square')
        
        if self.dimension == 3:
            ax = plt.axes(projection='3d')
            ax.scatter3D(coordinates[:, 1], coordinates[:, 0], coordinates[:, 2])
            print(self.LinkMatrix)
            ax.set_xlabel('x')
            ax.set_ylabel('t')
            ax.set_zlabel('y')

        plt.show()
            
    # CALCULATING LINK MATRIX
    
    def find_linkmatrix(self):
        #L_ij = 1 if e_i <* e_j; 0 otherwise, where <* is a link
        #L = C - f(C2)
        if self.LinkMatrix is None: 
            LinkMatrix: np.array = self.CausalMatrix - fC2(self.CausalMatrix)
            self.LinkMatrix = LinkMatrix
        
    
    # CALCULATING MYRHEIM_MEYER DIMENSION 
    
    def findOrderingFraction(self):
        '''finds ordering fraction of self.interval, r = 2R/ n(n-1)'''
        
        n = len(self.ElementList)
        r = 2*np.sum(self.CausalMatrix)/ (n*(n-1))
        return r 
    
    def Myhreim_Meyer_dimension(self, d, r): 
        #Solves d_MM(r) = 0
        
        if r < 0 or r > 1: 
            raise ValueError('Make sure ordering fraction, r is between 0 and 1!')
            
        return 1.5*gamma(d/2 + 1)*gamma(d + 1)/ gamma(3*d/2 + 1) - r
    
    def find_Myhreim_Meyer_dimension(self): 
        
        Myhreim_Meyer_dimension = fsolve(self.Myhreim_Meyer_dimension, x0 = 4, args = (self.findOrderingFraction()))
        # print(f'The Myhreim_Meyer_dimension is {Myhreim_Meyer_dimension}.')
        
        return Myhreim_Meyer_dimension
    
    # CALCULATING MOLECULES 
    
    def find_molecules(self):
        if self.LinkMatrix is None: 
            self.find_linkmatrix()
        maximals = []
        maximal_but_ones = []
        for i in range(len(self.LinkMatrix)):
            links = sum(self.LinkMatrix[i])
            if links == 0:
                maximals.append(i)
            elif links == 1:
                maximal_but_ones.append(i)
        
        H_array = []

        if not self.dynamic:
            for maximal in maximals:
                if self.ElementList[maximal].coordinates[0] > self.ElementList[maximal].coordinates[1]:
                    count = 0
                    for minimal_link in set(np.where(self.LinkMatrix[:, maximal] == 1)[0]).intersection(maximal_but_ones):
                        if self.ElementList[minimal_link].coordinates[0] < self.ElementList[minimal_link].coordinates[1]: 
                            count += 1 

                    while len(H_array) < count + 1:
                        H_array.append(0)
                    H_array[count] += 1
        else:
            for maximal in maximals:
                if inside_horizon(self.ElementList[maximal].coordinates):
                    count = 0
                    for minimal_link in set(np.where(self.LinkMatrix[:, maximal] == 1)[0]).intersection(maximal_but_ones):
                        if not inside_horizon(self.ElementList[minimal_link].coordinates): 
                            count += 1 

                    while len(H_array) < count + 1:
                        H_array.append(0)
                    H_array[count] += 1
        
        return H_array
    
if __name__ == "__main__":
       
    np.random.seed(0)

    tic = time.time()

    c = CausalSet(sprinkling_density = 0.1, 
                  dimension = 4, 
                  periodicBC = True,
                  DynamicBH = True,
                  T = 5)

    #c.visualisation()
    # print(c.ElementList)
    # print('Casual Matrix: \n', c.CausalMatrix)
    #c.find_linkmatrix()
    #print('Link Matrix: \n', c.LinkMatrix)
    #print('MM dimension is', c.find_Myhreim_Meyer_dimension())
    print('Number of Points:', len(c.ElementList))
    print(f'Spacetime Volume is {c.SpacetimeVolume}')
    #c.find_molecules()
    # c.visualisation()
    toc = time.time() 

    print(f'Time elapsed is {toc - tic}')
        
    
        
# =============================================================================
#     def main():     
#         np.random.seed(10)
# 
#         tic = time.time()
# 
#         c = CausalSet(sprinkling_density = 1000, 
#                       dimension = 3, 
#                       periodicBC = False)
#         # print(c.ElementList)
#         # print('Casual Matrix: \n', c.CausalMatrix)
#         #c.find_linkmatrix()
#         #print('Link Matrix: \n', c.LinkMatrix)
#         
#         print('MM dimension is', c.find_Myhreim_Meyer_dimension())
#         print('Number of Points:', len(c.ElementList))
#         #c.find_molecules()
#         # c.visualisation()
#         toc = time.time() 
#     
#         print(f'Time elapsed is {toc - tic}')
#         
#     #main()
# 
#     cProfile.run("main()", "output.dat")
#     
#     with open("output_time.txt", 'w') as f: 
#         p = pstats.Stats("output.dat", stream = f)
#         p.sort_stats("time").print_stats() 
#         
#     with open("output_calls.txt", "w") as f: 
#         p = pstats.Stats("output.dat", stream = f)
#         p.sort_stats("calls").print_stats()
# =============================================================================
        
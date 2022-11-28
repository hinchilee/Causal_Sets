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
from causalsetfunctions import spacetime_interval, inside_horizon
from causalEvent import CausalEvent
from Sprinkling import Sprinkling_Uniform, Sprinkling_Bicone

class CausalSet(object): 
    
    def __init__(self, **kwargs): 
        
        self.dimension = kwargs.get('dimension', 2)
        self.BHtype = kwargs.get('BHtype')
        self.ElementList: list (CausalEvent) = list() 
        
        #Sprinkling and sorting by time coordinate 
        
        # FOR RINDLER HORIZON  
        if self.BHtype == 'Rindler':
            # sprinkledcoords = Sprinkling_Bicone(dimension = self.dimension, number_of_points = kwargs.get('number_of_points', 5))
            # Normalised Sprinkling Volume to 1!!! Important to find out <N> then Poisson, not get N from Poisson then scale to area 
            
            self.SpacetimeVolume = 1
            sprinkledcoords = Sprinkling_Uniform(dimension = kwargs.get('dimension', 2),number_of_points = poisson.rvs(kwargs.get('sprinkling_density')), bounds = kwargs.get('bounds', np.array([[-0.5,0.5] for i in range(kwargs.get('dimension', 2))])))
            self.wrapAroundLength = 1

        # FOR DYNAMIC BLACKHOLE HORIZON
        elif self.BHtype == 'Dynamic':
            # Sprinkle uniformly into a slab to calculate entropy for Sigma_t = T 
            # into bounds of t: [T-1, T+1]; 
            # x, y, z = [-T-2, T+2]; wrap around y
            
            try: 
                self.T = kwargs.get('T')
            except: 
                raise ValueError('Please enter a value for T when sprinkling into a Dynamic BH!')
            spaceBounds = np.array([-self.T-2, self.T+2])
            timeBounds = np.array([self.T-1, self.T+1])
            self.SpacetimeVolume = ((2*(self.T+2))**(self.dimension - 1))*2
            AveragePoints = self.SpacetimeVolume*(kwargs.get('sprinkling_density'))
            boundsArray = np.concatenate((timeBounds, np.tile(spaceBounds, self.dimension - 1)),axis = 0).reshape(self.dimension, 2)
            sprinkledcoords = Sprinkling_Uniform(dimension = self.dimension, 
                                                number_of_points = poisson.rvs(AveragePoints),
                                                bounds = boundsArray) 
            self.wrapAroundLength = 2*(self.T+2)
        
        elif self.BHtype == 'Empty':
            
            self.SpacetimeVolume = 1
            sprinkledcoords = Sprinkling_Uniform(dimension = kwargs.get('dimension', 2),number_of_points = poisson.rvs(kwargs.get('sprinkling_density')), bounds = kwargs.get('bounds', np.array([[-0.5,0.5] for i in range(kwargs.get('dimension', 2))])))
            self.wrapAroundLength = 1
            
        else: 
            raise ValueError('Blackhole type should be Rindler, Dynamic or Empty!')
        
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
                    if spacetime_interval(self.ElementList[i].coordinates, self.ElementList[j].coordinates, self.BHtype, self.wrapAroundLength) < 0: 
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
            if self.LinkMatrix is None:
                self.find_linkmatrix()
            plt.figure(figsize= (12,8))
            plt.scatter(coordinates[:, 1], coordinates[:, 0], s = 70, c = 'black')
            for i in range(len(self.LinkMatrix)):
                for j in range(len(self.LinkMatrix[i])):
                    if self.LinkMatrix[i][j] == 1:
                        plt.plot([coordinates[i][1], coordinates[j][1]], [coordinates[i][0], coordinates[j][0]], color = 'royalblue',linewidth = 0.8)
            plt.xlabel('Space', fontsize = 20 )
            plt.ylabel('Time', fontsize = 20 )
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
            LinkMatrix: np.array = self.CausalMatrix - np.matmul(self.CausalMatrix, self.CausalMatrix).clip(0, 1)
            
            # LinkMatrix = np.zeros(self.CausalMatrix.shape, dtype = int)
            # n = len(self.ElementList)
            
            # for i in range(n): 
            #     for j in range(i+1, n): 
            #         # 1. Check Cij = 1 
            #         if self.CausalMatrix[i][j] == 1: 
            #             # 2. Check C^2_ij = 0 
            #             if np.sum(self.CausalMatrix[i,:]*self.CausalMatrix[:, j]) == 0: 
            #                 LinkMatrix[i][j] = 1
                    
            #             else:
            #                 LinkMatrix[i][j] = 0
            #         else: 
            #             LinkMatrix[i][j] = 0       
            
            # LinkMatrix = np.zeros(self.CausalMatrix.shape, dtype = int)
            # n = len(self.ElementList)
            
            # for i in range(n): 
            #     for j in range(i+1, n): 
            #         # 1. Check Cij = 1 
            #         if self.CausalMatrix[i][j] == 1: 
            #             # 2. Check C^2_ij = 0 
            #             for k in range(n):
            #                   if self.CausalMatrix[i][k]*self.CausalMatrix[k][j] == 0: 
            #                       pass 
            #                   else:
            #                       break
            #                   LinkMatrix[i][j] = 1
            #         else: 
            #             LinkMatrix[i][j] = 0        
            
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

        if self.BHtype == 'Rindler':
            for maximal in maximals:
                if self.ElementList[maximal].coordinates[0] > self.ElementList[maximal].coordinates[1]:
                    count = 0
                    for minimal_link in set(np.where(self.LinkMatrix[:, maximal] == 1)[0]).intersection(maximal_but_ones):
                        if self.ElementList[minimal_link].coordinates[0] < self.ElementList[minimal_link].coordinates[1]: 
                            count += 1 

                    while len(H_array) < count + 1:
                        H_array.append(0)
                    H_array[count] += 1
        elif self.BHtype == 'Dynamic':
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
      

    def main():     
        np.random.seed(0)

        tic = time.time()
    
        c = CausalSet(sprinkling_density = 1000,    # 0.1-1 for Dynamic, 1k - 10k for Rindler, Empty 
                      dimension = 4, 
                      BHtype = 'Rindler',           # 'Rindler', 'Dynamic', 'Empty' 
                      T = 5)                        # T is only needed when BHtype = 'Dynamic'
    
        #c.visualisation()
        # print(c.ElementList)
        # print('Casual Matrix: \n', c.CausalMatrix)
        # c.find_linkmatrix()
        #print('Link Matrix: \n', c.LinkMatrix)
        #print('MM dimension is', c.find_Myhreim_Meyer_dimension())
        print('Number of Points:', len(c.ElementList))
        print(f'Spacetime Volume is {c.SpacetimeVolume}')
        print(c.find_molecules())
        #c.visualisation()
        toc = time.time() 
    
        print(f'Time elapsed is {toc - tic}')
    

    cProfile.run("main()", "output.dat")
    
    with open("output_time.txt", 'w') as f: 
        p = pstats.Stats("output.dat", stream = f)
        p.sort_stats("time").print_stats() 
        
    with open("output_calls.txt", "w") as f: 
        p = pstats.Stats("output.dat", stream = f)
        p.sort_stats("calls").print_stats()
        
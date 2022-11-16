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
from itertools import combinations 
from causalsetfunctions import f, findC2, spacetime_interval, generate_CausalMatrix
from causalEvent import CausalEvent
from Sprinkling import Sprinkling_Uniform, Sprinkling_Bicone
import multiprocessing as mp

class CausalSet(object): 
    
    def __init__(self, **kwargs): 
        
        self.dimension = kwargs.get('dimension', 2)
        self.ElementList: list (CausalEvent) = list() 
    
        #Sprinkling and sorting by time coordinate 
        
# =============================================================================
#         #2d uniform distribution
#         sprinkledcoords = Sprinkling_Uniform(dimension = kwargs.get('dimension', 2), number_of_points = kwargs.get('number_of_points', 5), bounds = kwargs.get('bounds', np.array([[0,1] for i in range(kwargs.get('dimension', 2))])))
#         #Manually add two top and bottom points
#         sprinkledcoords = np.concatenate((sprinkledcoords, np.array([0,0]).reshape(1,2)))
#         sprinkledcoords = np.concatenate((sprinkledcoords, np.array([1,1]).reshape(1,2)))
#         t = (sprinkledcoords[:, 1] + sprinkledcoords[:, 0])/np.sqrt(2)
#         x = (sprinkledcoords[:, 1] - sprinkledcoords[:, 0])/np.sqrt(2)
#         sprinkledcoords[:, 0] = t
#         sprinkledcoords[:, 1] = x
# =============================================================================
        
        # sprinkledcoords = Sprinkling_Bicone(dimension = self.dimension, number_of_points = kwargs.get('number_of_points', 5))
        sprinkledcoords = Sprinkling_Uniform(dimension = kwargs.get('dimension', 2), number_of_points = poisson.rvs(kwargs.get('sprinkling_density', 5)), bounds = kwargs.get('bounds', np.array([[-1,1] for i in range(kwargs.get('dimension', 2))])))

        #sort by time
        sprinkledcoords = sprinkledcoords[sprinkledcoords[:, 0].argsort()]
        
        for i, coords in enumerate(sprinkledcoords): 
            self.ElementList.append(CausalEvent(label = i, coordinates = coords))
    
        #Create causal matrix using spacetime intervals
        
# =============================================================================
#         tic3 = time.time()
#         self.CausalMatrix: np.array = np.zeros((len(self.ElementList), len(self.ElementList)),dtype = int)
#         for i, CausalEvent1 in enumerate(self.ElementList):
#             for j, CausalEvent2 in enumerate(self.ElementList[i:]):
#                 #if CausalEvent1.label < CausalEvent2.label: 
#                 self.CausalMatrix[i, j+i] = 1 if spacetime_interval(CausalEvent1.coordinates, CausalEvent2.coordinates) < 0 else 0
#         toc3 = time.time()
#         print(f'time taken to generate causal matrix by pairwise is {toc3 - tic3}')
#     
# =============================================================================
        self.LinkMatrix = None
        self.Interval: set(CausalEvent) = set()
        tic2 = time.time() 
        
        ElementList = [[ele.coordinates for ele in self.ElementList]]

        # Intialise Parallel Processing
        # pool = mp.Pool(mp.cpu_count() - 4) # leave 4 cores unused to avoid overloading
        # self.CausalMatrix = pool.map(generate_CausalMatrix, ElementList)[0]
        # pool.close() 
        
        self.CausalMatrix = self.generate_CausalMatrix()
        toc2 = time.time() 
        # print(f'time taken to generate causal matrix by transitivity is {toc2 - tic2}')
        
    def generate_CausalMatrix(self): 
        # Induce causal relations by transitivity
        
        A = np.zeros((len(self.ElementList), len(self.ElementList)), dtype = int)
        for j in range(len(self.ElementList)): 
            for i in reversed(range(j)): 
                if A[i,j] == 0:
                    if spacetime_interval(self.ElementList[j].coordinates, self.ElementList[i].coordinates) < 0: 
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
            
    #CALCULATING LINK MATRIX
    
    def find_linkmatrix(self):
        #L_ij = 1 if e_i <* e_j; 0 otherwise, where <* is a link
        #L = C - f(C2)
        if self.LinkMatrix is None: 
            LinkMatrix: np.array = self.CausalMatrix - f(findC2(self.CausalMatrix))
            self.LinkMatrix = LinkMatrix
        
        return self.LinkMatrix
    
    #CALCULATING MYRHEIM_MEYER DIMENSION 
    
    def find_interval(self, CausalEvent1, CausalEvent2):
        '''Calculates interval given two causal events and updates it to self.Interval'''
        
        if self.to_the_past(CausalEvent1, CausalEvent2): 
            CEpast, CEfuture = CausalEvent1.label, CausalEvent2.label
        elif self.to_the_future(CausalEvent1, CausalEvent2): 
            CEfuture, CEpast = CausalEvent1.label, CausalEvent2.label
        
        interval: set(CausalEvent) = CEpast.future & CEfuture.past
        self.Interval = interval
        
        return interval
    
    def findOrderingFraction(self):
        '''finds ordering fraction of self.interval, r = 2R/ n(n-1)'''
        
# =============================================================================
#         #get element labels 
#         #intervalLabels = [CausalEvent.label for CausalEvent in self.interval]
#         intervalLabels = [CausalEvent.label for CausalEvent in self.ElementList]
#         #generate a pairwise list
#         pairs = list(combinations(intervalLabels, 2))
#         RelationsCount = 0 
#         for pair in pairs: 
#             RelationsCount += self.CausalMatrix[pair[0], pair[1]]
#             RelationsCount += self.CausalMatrix[pair[1], pair[0]]
#         
#         #n = len(self.interval)
#         
#         
# =============================================================================
        
    
        
        n = len(self.ElementList)
        #r = 2*RelationsCount/ (n*(n-1))
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
    
    #CHECKING ELEMENTAL RELATIONS 
    
    def to_the_past(self, CausalEvent1, CausalEvent2): 
        '''Checks whether ele1 < ele2'''
    
        l1 = CausalEvent1.label 
        l2 = CausalEvent2.label 
        
        if self.CausalMatrix[l1][l2] == 1:
            return True
        else: 
            return False
    
    def to_the_future(self, CausalEvent1, CausalEvent2): 
        '''Checks whether ele1 > ele2'''
        
        l1 = CausalEvent1.label 
        l2 = CausalEvent2.label 
        
        if self.CausalMatrix[l2][l1] == 1:
            return True
        else: 
            return False
        
    def spacelike_to(self, CausalEvent1, CausalEvent2): 
        '''Checks whether ele1 is spacelike to ele2'''
        
        l1 = CausalEvent1.label 
        l2 = CausalEvent2.label 
        
        if self.CausalMatrix[l1][l2] == 0 and self.CausalMatrix[l2][l1] == 0:
            return True
        else: 
            return False
        
        
    def check_relation(self, CausalEvent1, CausalEvent2): 
        '''gives relation R in ele1 R ele2, where R = {precedes, succeeds, spacelike to}'''
        
        l1 = CausalEvent1.label 
        l2 = CausalEvent2.label 
        if self.CausalMatrix[l1][l2] == 1:
            state = 'precedes'
        elif self.CausalMatrix[l2][l1] == 1: 
            state = 'succedes'
        else: 
            state = 'is spacelike to'
        print(f'Causal Event {l1} {state} Causal Event {l2}.')

    def find_molecules(self):
        self.find_linkmatrix()
        maximals = []
        maximal_but_ones = []
        two_chains = []
        for i in range(len(self.LinkMatrix)):
            links = sum(self.LinkMatrix[i])
            if links == 0:
                maximals.append(i)
            elif links == 1:
                maximal_but_ones.append(i)
        
        for maximal in maximals:
            for minimal_link in set(np.where(self.LinkMatrix[:, maximal] == 1)[0]).intersection(maximal_but_ones):
                two_chains.append((minimal_link, maximal))

        count = 0
        for chain in two_chains:
            if self.ElementList[chain[0]].coordinates[0] < self.ElementList[chain[0]].coordinates[1] and self.ElementList[chain[1]].coordinates[0] > self.ElementList[chain[1]].coordinates[1]:
                count += 1

        return count
        
if __name__ == "__main__":
    
    def main():     
        np.random.seed(10)

        tic = time.time()

        c = CausalSet(sprinkling_density = 1000, dimension = 4)
        # print(c.ElementList)
        # print('Casual Matrix: \n', c.CausalMatrix)
        # c.find_linkmatrix()
        # print('Link Matrix: \n', c.LinkMatrix)
        
        # print(c.find_Myhreim_Meyer_dimension())
        print('Number of Points:', len(c.ElementList))
        c.find_molecules()
        # c.visualisation()
        toc = time.time() 
    
        print(f'Time elapsed is {toc - tic}')
        
    #main()

    cProfile.run("main()", "output.dat")
    
    with open("output_time.txt", 'w') as f: 
        p = pstats.Stats("output.dat", stream = f)
        p.sort_stats("time").print_stats() 
        
    with open("output_calls.txt", "w") as f: 
        p = pstats.Stats("output.dat", stream = f)
        p.sort_stats("calls").print_stats()
        
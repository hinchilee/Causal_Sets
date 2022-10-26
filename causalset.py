# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 12:09:46 2022

@author: leehi
"""


import numpy as np 
import time 
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import gamma 
from itertools import combinations 
from causalsetfunctions import f, findC2, spacetime_interval
from causalEvent import CausalEvent
from Sprinkling import Sprinkling_Minkowski 

np.random.seed(10)

class CausalSet(object): 
    
    def __init__(self, **kwargs): 
        
        self.ElementList: list (CausalEvent) = list() 
        #Sprinkling and sorting by time coordinate 
        sprinkledcoords = Sprinkling_Minkowski(dimension = kwargs.get('dimension', 2), number_of_points = kwargs.get('number_of_points', 5), bounds = kwargs.get('bounds', np.array([[0,1] for i in range(kwargs.get('dimension', 2))])))
        t = (sprinkledcoords[:, 1] + sprinkledcoords[:, 0])/np.sqrt(2)
        x = (sprinkledcoords[:, 1] - sprinkledcoords[:, 0])/np.sqrt(2)
        sprinkledcoords[:, 0] = t
        sprinkledcoords[:, 1] = x
        sprinkledcoords = sprinkledcoords[sprinkledcoords[:, 0].argsort()]
        for i, coords in enumerate(sprinkledcoords): 
            self.ElementList.append(CausalEvent(label = i, coordinates = coords))
    
        #Create causal matrix using spacetime intervals
        self.CausalMatrix: np.array = np.zeros((len(self.ElementList), len(self.ElementList)))
        for CausalEvent1 in self.ElementList:
            for CausalEvent2 in self.ElementList:
                self.CausalMatrix[CausalEvent1.label, CausalEvent2.label] = 1 if (spacetime_interval(CausalEvent1.coordinates, CausalEvent2.coordinates) < 0 and CausalEvent1.label < CausalEvent2.label) else 0
        
        self.LinkMatrix: np.array = self.find_linkmatrix()
        self.Interval: set(CausalEvent) = set() 
        
        #TODO Add future and past of causal events
        
    #VISUALIATION 
    
    def visualisation(self): 
    
        coordinates = np.array([x.coordinates for x in self.ElementList])
        plt.scatter(coordinates[:, 1], coordinates[:, 0], s = 100)
        # U = (coordinates[:, 0] + coordinates[:, 1])/2
        # V = (coordinates[:, 0] - coordinates[:, 1])/2
        # plt.scatter(U, V, s=100)

        for i in range(len(self.LinkMatrix)):
            for j in range(len(self.LinkMatrix[i])):
                if self.LinkMatrix[i][j] == 1:
                    # plt.plot([U[i], U[j]], [V[i], V[j]], color = 'green')
                    plt.plot([coordinates[i][1], coordinates[j][1]], [coordinates[i][0], coordinates[j][0]], color = 'green')
        plt.xlabel('Space', fontsize = 40)
        plt.ylabel('Time', fontsize = 40)
        plt.xticks(fontsize = 30)
        plt.yticks(fontsize = 30)
        plt.show()
            
    #CALCULATING LINK MATRIX
    
    def find_linkmatrix(self):
        #L_ij = 1 if e_i <* e_j; 0 otherwise, where <* is a link
        #L = C - f(C2)
        
        LinkMatrix: np.array = self.CausalMatrix - f(findC2(self.CausalMatrix))
        
        return LinkMatrix
    
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
        
        #get element labels 
        intervalLabels = [CausalEvent.label for CausalEvent in self.interval]
        #generate a pairwise list
        pairs = list(combinations(intervalLabels, 2))
        RelationsCount = 0 
        for pair in pairs: 
            RelationsCount += self.CausalMatrix[pair[0], pair[1]]
            RelationsCount += self.CausalMatrix[pair[1], pair[0]]
        
        n = len(self.interval)
        r = 2*RelationsCount/ (n*(n-1))
        return r 
    
    def Myhreim_Meyer_dimension(self, d, r): 
        #Solves d_MM(r) = 0
        
        if r < 0 or r > 1: 
            raise ValueError('Make sure ordering fraction, r is between 0 and 1!')
            
        return 1.5*gamma(d/2)*gamma(d)/ gamma(3*d/2) - r
    
    def find_Myhreim_Meyer_dimension(self): 
        
        Myhreim_Meyer_dimension = fsolve(self.Myhreim_Meyer_dimension, x0 = 4, args = (self.findOrderingFraction))
        print(f'The Myhreim_Meyer_dimension is {Myhreim_Meyer_dimension}.')
        
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
        
if __name__ == "__main__":
    tic = time.time()
    c = CausalSet(number_of_points = 100) 
    # print(c.ElementList)
    # print('Casual Matrix: \n', c.CausalMatrix)
    # print('Link Matrix: \n', c.LinkMatrix)
    c.visualisation()
    toc = time.time() 
    print(f'Time elapsed is {toc - tic}')
        
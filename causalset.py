# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 12:09:46 2022

@author: leehi
"""

import cProfile
import pstats
import io
from pstats import SortKey
import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.optimize import fsolve
from scipy.special import gamma
from causalsetfunctions import spacetime_interval, inside_horizon, n_ball_volume
from causalEvent import CausalEvent
from Sprinkling import Sprinkling_Uniform, Sprinkling_Bicone, Sprinkling_Tube
import multiprocessing as mp


class CausalSet(object):

    def __init__(self, **kwargs):

        self.dimension = kwargs.get('dimension', 2)
        self.BHtype = kwargs.get('BHtype')
        self.sprinkling = kwargs.get('sprinkling', 'Uniform')
        self.ElementList: list(CausalEvent) = list()

        # Sprinkling and sorting by time coordinate

        # FOR RINDLER HORIZON
        if self.BHtype == 'Rindler':
            # sprinkledcoords = Sprinkling_Bicone(dimension = self.dimension, number_of_points = kwargs.get('number_of_points', 5))
            # Normalised Sprinkling Volume to 1!!! Important to find out <N> then Poisson, not get N from Poisson then scale to area

            bounds = kwargs.get('bounds', np.array([[-0.5, 0.5] for i in range(kwargs.get('dimension', 2))]))
            #print('Boundaries', bounds)
            self.SpacetimeVolume = np.prod(bounds[:,1] - bounds[:,0])
            print('Spacetime volume', self.SpacetimeVolume)
            AveragePoints = self.SpacetimeVolume * (kwargs.get('sprinkling_density'))
            noPoints = poisson.rvs(AveragePoints)
            print('Number of points:', noPoints)
            sprinkledcoords = Sprinkling_Uniform(dimension=kwargs.get('dimension', 2), number_of_points=noPoints, bounds=bounds)
            self.wrapAroundLength = 1

        # FOR DYNAMIC BLACKHOLE HORIZON
        elif self.BHtype == 'Dynamic':
            # Sprinkle uniformly into a slab to calculate entropy for Sigma_t = T
            # into bounds of t: [T-1, T+1];
            # x, y, z = [-T-2, T+2]; wrap around y
            self.wrapAroundLength = 1  #Placeholder - isn't used
            self.T = kwargs.get('T', 3)
            
            if self.sprinkling == 'Uniform':
                spaceBounds = np.array([-self.T-2, self.T+2])
                timeBounds = np.array([self.T-3, self.T])
                bounds = np.concatenate((timeBounds, np.tile(
                    spaceBounds, self.dimension - 1)), axis=0).reshape(self.dimension, 2)
                self.SpacetimeVolume = np.prod(bounds[:,1] - bounds[:,0])
                print('Spacetime volume', self.SpacetimeVolume)
                AveragePoints = self.SpacetimeVolume * (kwargs.get('sprinkling_density'))
                noPoints = poisson.rvs(AveragePoints)
                print('Number of points:', noPoints)
                sprinkledcoords = Sprinkling_Uniform(dimension=self.dimension,
                                                     number_of_points=poisson.rvs(
                                                         noPoints),
                                                     bounds=bounds)
                
            elif self.sprinkling == 'Tube':
                ndimension = self.dimension - 1
                try:
                    R_min, R_max, T_min, T_max = kwargs.get('bounds', [self.T * 0.7, self.T + 1, self.T - 3, self.T]) 
                except: 
                    raise ValueError('Correct bounds should be in the form of [R_min, R_max, T_min, T_max]')
                self.SpacetimeVolume = (T_max-T_min)*(n_ball_volume(ndimension, R_max) - n_ball_volume(ndimension, R_min))
                print('Spacetime volume', self.SpacetimeVolume)
                AveragePoints = self.SpacetimeVolume * (kwargs.get('sprinkling_density'))
                noPoints = poisson.rvs(AveragePoints)
                print('Number of points:', noPoints)
                sprinkledcoords = Sprinkling_Tube(dimension=self.dimension,
                                                  number_of_points=noPoints,
                                                  bounds = [R_min, R_max, T_min, T_max])

        elif self.BHtype == 'Empty':

            self.SpacetimeVolume = 1
            # sprinkledcoords = Sprinkling_Uniform(dimension = kwargs.get('dimension', 2),number_of_points = poisson.rvs(kwargs.get('sprinkling_density')), bounds = kwargs.get('bounds', np.array([[-0.5,0.5] for i in range(kwargs.get('dimension', 2))])))
            # sprinkledcoords = Sprinkling_Bicone(dimension = self.dimension, number_of_points = poisson.rvs(kwargs.get('sprinkling_density')))
            sprinkledcoords = Sprinkling_Tube(dimension=self.dimension, number_of_points=poisson.rvs(
                kwargs.get('sprinkling_density')), R_min=0.5, R_max=1, T_min=0, T_max=1)
            self.wrapAroundLength = 1

        else:
            raise ValueError(
                'Blackhole type should be Rindler, Dynamic or Empty!')

        # sort by time
        sprinkledcoords = sprinkledcoords[sprinkledcoords[:, 0].argsort()]

        for i, coords in enumerate(sprinkledcoords):
            self.ElementList.append(CausalEvent(label=i, coordinates=coords))

        self.LinkMatrix = None
        self.VElementsLabelsList = list()
        self.generate_CausalMatrix()

    def generate_CausalMatrix(self):
        # Induce causal relations by transitivity

        # A = np.zeros((len(self.ElementList), len(self.ElementList)), dtype = int)
        # for j in range(len(self.ElementList)):
        #     for i in reversed(range(j)):
        #         if A[i,j] == 0:
        #             if spacetime_interval(self.ElementList[i].coordinates, self.ElementList[j].coordinates, self.BHtype, self.wrapAroundLength) < 0:
        #                 A[i,j] = 1
        #                 #Then inherit i's past
        #                 A[:,j] = np.bitwise_or(A[:,j], A[:,i])
        #             else:
        #                 pass
        #         else:
        #             pass
        # self.CausalMatrix = A

        self.cores = mp.cpu_count() - 3
        # print(f'Cores: {self.cores}')
        self.CausalMatrix = np.zeros(
            (len(self.ElementList), len(self.ElementList)), dtype=int)
        with mp.Pool(self.cores) as pool:
            for resultrow, result in enumerate(pool.map(self.causalmatrixpooltask, range(len(self.ElementList)))):
                self.CausalMatrix[resultrow, (resultrow+1):] = result

    def causalmatrixpooltask(self, i):
        res = list()
        for j in range(i+1, len(self.ElementList)):
            if spacetime_interval(self.ElementList[i].coordinates, self.ElementList[j].coordinates, self.BHtype, self.wrapAroundLength) < 0:
                res.append(1)
            else:
                res.append(0)
        return res

    # VISUALIATION

    def visualisation(self):
        coordinates = np.array([x.coordinates for x in self.ElementList])

        if self.dimension == 2:
            if self.LinkMatrix is None:
                self.find_linkmatrix()
            plt.figure(figsize=(12, 8))
            plt.scatter(coordinates[:, 1], coordinates[:, 0], s=70, c='black')
            for i in range(len(self.LinkMatrix)):
                for j in range(len(self.LinkMatrix[i])):
                    if self.LinkMatrix[i][j] == 1:
                        plt.plot([coordinates[i][1], coordinates[j][1]], [
                                 coordinates[i][0], coordinates[j][0]], color='royalblue', linewidth=0.8)
            if self.BHtype == 'Rindler':
                xlinspace = np.linspace(-0.5, 0.5, 100)
                plt.plot(xlinspace, xlinspace,
                         label='Rindler Horizon', c='red')
            elif self.BHtype == 'Dynamic':
                xlinspace = np.linspace(0, 10, 100)
                plt.plot(xlinspace, xlinspace,
                         label='Dynamic Horizon', c='red')
                plt.plot(-xlinspace, xlinspace,
                         label='Dynamic Horizon', c='red')
                xlinspace2 = np.linspace(-10, 10, 200)
                plt.plot(xlinspace2, [self.T]*len(xlinspace2),
                         label=f'Sigma plane t = {self.T}', c='green')
            plt.xlabel('Space', fontsize=20)
            plt.ylabel('Time', fontsize=20)
            plt.axis('square')
            plt.legend()

        if self.dimension == 3:
            ax = plt.axes(projection='3d')
            # ?, t, ?
            ax.scatter3D(coordinates[:, 1],
                         coordinates[:, 0], coordinates[:, 2])
            print(self.LinkMatrix)
            ax.set_xlabel('x')
            ax.set_ylabel('t')
            ax.set_zlabel('y')

        plt.show()

    # CALCULATING LINK MATRIX

    def find_linkmatrix(self):
        # L_ij = 1 if e_i <* e_j; 0 otherwise, where <* is a link
        # L = C - f(C2)
        if self.LinkMatrix is None:
            # LinkMatrix: np.array = self.CausalMatrix - np.matmul(self.CausalMatrix, self.CausalMatrix).clip(0, 1)

            self.LinkMatrix = np.zeros(self.CausalMatrix.shape, dtype=int)
            with mp.Pool(self.cores) as pool:
                for resultrow, result in enumerate(pool.map(self.linkmatrixpooltask, range(len(self.ElementList)))):
                    self.LinkMatrix[resultrow, (resultrow+1):] = result

    def linkmatrixpooltask(self, i):
        res = list()
        for j in range(i+1, len(self.ElementList)):
            # 1. Check Cij = 1
            if self.CausalMatrix[i][j] == 1:
                # 2. Check C^2_ij = 0
                if np.sum(self.CausalMatrix[i, :]*self.CausalMatrix[:, j]) == 0:
                    res.append(1)
                else:
                    res.append(0)
            else:
                res.append(0)
        return res

    # CALCULATING MYRHEIM_MEYER DIMENSION

    def findOrderingFraction(self):
        '''finds ordering fraction of self.interval, r = 2R/ n(n-1)'''

        n = len(self.ElementList)
        r = 2*np.sum(self.CausalMatrix) / (n*(n-1))
        return r

    def Myhreim_Meyer_dimension(self, d, r):
        # Solves d_MM(r) = 0

        if r < 0 or r > 1:
            raise ValueError(
                'Make sure ordering fraction, r is between 0 and 1!')

        return 1.5*gamma(d/2 + 1)*gamma(d + 1) / gamma(3*d/2 + 1) - r

    def find_Myhreim_Meyer_dimension(self):

        Myhreim_Meyer_dimension = fsolve(
            self.Myhreim_Meyer_dimension, x0=4, args=(self.findOrderingFraction()))
        # print(f'The Myhreim_Meyer_dimension is {Myhreim_Meyer_dimension}.')

        return Myhreim_Meyer_dimension

    # CALCULATING MOLECULES

    def find_molecules(self):
        # if self.LinkMatrix is None:
        #     self.find_linkmatrix()
        maximals = []
        maximal_but_ones = []

        if self.BHtype == 'Rindler':
            for i in range(len(self.CausalMatrix)):
                links = sum(self.CausalMatrix[i])
                if links == 0:
                    maximals.append(i)
                elif links == 1:
                    maximal_but_ones.append(i)

        elif self.BHtype == 'Dynamic':
            diffarray = np.array([ele.coordinates[0]
                                 for ele in self.ElementList]) - self.T
            try:
                cutoffindex = np.min(np.where(diffarray > 0))
                # Crop Causal matrix so that only elements before sigmaT are included
                self.croppedCausalMatrix = self.CausalMatrix[:cutoffindex, :cutoffindex]

            except ValueError:
                self.croppedCausalMatrix = self.CausalMatrix

            for i in range(len(self.croppedCausalMatrix)):
                links = sum(self.croppedCausalMatrix[i])
                if links == 0:
                    maximals.append(i)
                elif links == 1:
                    maximal_but_ones.append(i)

        H_array = []
        min_time = 10000
        min_distance = 10000 
        max_distance = -10000
        if self.BHtype == 'Rindler':
    
            for maximal in maximals:
                if self.ElementList[maximal].coordinates[0] > self.ElementList[maximal].coordinates[1]:
                    count = 0
                    for minimal_link in set(np.where(self.CausalMatrix[:, maximal] == 1)[0]).intersection(maximal_but_ones):
                        if self.ElementList[minimal_link].coordinates[0] < self.ElementList[minimal_link].coordinates[1]:
                            count += 1
                            min_time = min(min_time, (self.ElementList[minimal_link].coordinates[0] - 0.5))
                            min_distance = min([min_distance, self.ElementList[minimal_link].coordinates[1] - abs(self.ElementList[minimal_link].coordinates[0])+ 0.5])
                            max_distance = max([max_distance, self.ElementList[minimal_link].coordinates[1] + abs(self.ElementList[minimal_link].coordinates[0])+ 0.5])

                    while len(H_array) < count:
                        H_array.append(0)
                    if count > 0:
                        H_array[count - 1] += 1

        elif self.BHtype == 'Dynamic':
            for maximal in maximals:
                if inside_horizon(self.ElementList[maximal].coordinates):
                    count = 0
                    for minimal_link in set(np.where(self.croppedCausalMatrix[:, maximal] == 1)[0]).intersection(maximal_but_ones):
                        if not inside_horizon(self.ElementList[minimal_link].coordinates):
                            count += 1
                            min_time = min(
                                min_time, (self.ElementList[minimal_link].coordinates[0] - self.T))
                            minimallinknorm = np.linalg.norm(self.ElementList[minimal_link].coordinates[1:])
                            #Drag out light cone from maximal but one
                            min_distance = min([min_distance, (minimallinknorm - abs(self.ElementList[minimal_link].coordinates[0]))])
                            max_distance = max([max_distance, (minimallinknorm + abs(self.ElementList[minimal_link].coordinates[0]))])

                    while len(H_array) < count:
                        H_array.append(0)
                    if count > 0:
                        H_array[count - 1] += 1
        try:
            self.min_time = min_time
            self.max_distance = max_distance 
            self.min_distance = min_distance    #for dynamic, min distance is min spatial distance from the t-axis
            print(f'The minimum time of a molecule is {min_time}')
            print(f'The maximum distance of a molecule is {max_distance}')
            print(f'The minimum distance of a molecule is {min_distance}')
        except: 
            pass

        return H_array
    
    def find_Vmolecules(self):
        
        VElementsLabelsSet = set()
        V_count = 0 
        
        # 1. Check lower element has 2 future elements 
        for i in range(len(self.CausalMatrix)):
            row = self.CausalMatrix[i]
            linksBottomElement = sum(row)
            inside = False
            outside = False
            
            if linksBottomElement == 2:
                
                future2Elements = np.nonzero(row)
                
                # 2. Check both top elements has no causal future
                for j in future2Elements[0]: 
                    linksTopElements = sum(self.CausalMatrix[j])
                    if linksTopElements != 0: 
                        break 
                    
                    # 3. Check each element is outside H and inside H respectively 
                    if self.BHtype == 'Rindler': 
                        if self.ElementList[j].coordinates[0] > self.ElementList[j].coordinates[1]: 
                            inside = True
                        elif self.ElementList[j].coordinates[1] > self.ElementList[j].coordinates[0]:
                            outside = True 
                    elif self.BHtype == 'Dynamic': 
                        if inside_horizon(self.ElementList[j].coordinates): 
                            inside = True
                        elif inside_horizon(self.ElementList[j].coordinates) == False:
                            outside = True
                          
            # 4. Succuesfully identify V-molecule
            if (inside, outside) == (True, True): 
                V_count += 1 
                VElementsLabelsSet.add(future2Elements[0][0]) 
                VElementsLabelsSet.add(future2Elements[0][1]) 
                VElementsLabelsSet.add(i)
                
        self.VElementsLabelsList = list(VElementsLabelsSet)
        
        print(f'V-molecule count: {V_count}')
        
        return V_count
    
    


if __name__ == "__main__":


    #def main():
    np.random.seed(12)

    tic = time.time()
    boundsArray = np.array([[-0.5, 0.5] for i in range(3)])
    boundsArray[0][1] = 0.5 #no normalisation
    boundsArray[1][1] = 1.5

    c = CausalSet(sprinkling_density=500,    # 0.1-1 for Dynamic Uniform, 1k - 10k for Dynamic Tube, 1k - 10k for Rindler, Empty
                  dimension=4,
                  BHtype='Rindler',           # 'Rindler', 'Dynamic', 'Empty'
                  sprinkling='Tube',     # 'Uniform' or 'Tube' for 'Dynamic'BH 
                  T = 3)#,                    # T is only needed when BHtype = 'Dynamic'
                  #bounds = [5, 10, 2, 3])   # bounds for tube sprinkling in the form of [R_min, R_max, T_min, T_max]            

    #c.visualisation()
    # print(c.ElementList)
    #print('Casual Matrix: \n', c.CausalMatrix)
    # C2 = c.CausalMatrix
    # c.find_linkmatrix()
    # print('MM dimension is', c.find_Myhreim_Meyer_dimension())

    #print(c.find_molecules())
    c.find_Vmolecules()
    # print('Link Matrix: \n', c.LinkMatrix)
    #c.visualisation()

    toc = time.time()

    print(f'Time elapsed is {toc - tic}')
    
    
    # cProfile.run("main()", "output.dat")

    # with open("output_time.txt", 'w') as f:
    #     p = pstats.Stats("output.dat", stream = f)
    #     p.sort_stats("time").print_stats()

    # with open("output_calls.txt", "w") as f:
    #     p = pstats.Stats("output.dat", stream = f)
    #     p.sort_stats("calls").print_stats()

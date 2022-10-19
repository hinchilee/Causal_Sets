
import numpy as np 
from scipy.optimize import fsolve
from scipy.special import gamma 
from itertools import combinations 
from causalEvent import CausalEvent

class CausalSet(object): 
    
    def __init__(self, **kwargs): 
        #Sprinkling to be done 
        #Labelling elements 
        #Creation of Causal Sets
        
        self.Set: set(CausalEvent) = kwargs.get("Set")
        self.CausalMatrix: np.array = kwargs.get("CausalMatrix")
        self.LinkMatrix: np.array = kwargs.get("LinkMatrix") 
        self.Interval: set(CausalEvent) = set() 
        
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
        
            
        
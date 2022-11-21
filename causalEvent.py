
import numpy as np 

class CausalEvent(object):
    
    def __init__(self, label, **kwargs):
        self.coordinates: np.ndarray = kwargs.get('coordinates')
        self.label = label 
        self.insideHorizon = kwargs.get('insideHorizon')
    
    def __repr__(self): 
        return f'Causal event {self.label} at coordinates {self.coordinates}'
        
        
        
        
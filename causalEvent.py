
import numpy as np 

class CausalEvent(object):
    
    def __init__(self, label, **kwargs):
        self.past: set() = kwargs.get('past')
        self.future: set() = kwargs.get('future') 
        self.coordinates: np.ndarray = kwargs.get('coordinates')
        #self.label: int = kwargs.get('label') 
        self.label = label 
        
        
        
        
        

import numpy as np 

class CausalEvent(object):
    
    def __init__(self, label, **kwargs):
        # self.past = set()
        # self.future = set()
        self.coordinates: np.ndarray = kwargs.get('coordinates')
        #self.label: int = kwargs.get('label') 
        self.label = label 
    
    def __repr__(self): 
        return f'Causal event {self.label} at coordinates {self.coordinates}'
        
        
        
        
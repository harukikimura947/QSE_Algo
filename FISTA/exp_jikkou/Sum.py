import numpy as np
import add_king

class Sum:
    
    def __init__(self, star, dark):
        
        self.star = star
        self.dark = dark
        
        
    def add_1(self, a, b):
        
        return a + b

    def add_2(self, c, d):
        
        return c + d
    
    def add(self, a, b, c, d):
        
        result = add_king.add_king(self.add_1, self.add_2, a, b, c, d)
        
        return result
    
    
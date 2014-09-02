"""
Created on 6/8/2014

@author: victor
"""
from pyproct.data.matrix.matrixCalculator import MatrixCalculator

class DataDriver(object):
    
    def __init__(self, parameters):
        self.data_type = parameters["type"]
        self.data_handler = self.__create_data_handler(parameters)
        self.matrix_handler = MatrixCalculator.calculate(self.data_handler, 
                                                         parameters["matrix"])
        
        
    
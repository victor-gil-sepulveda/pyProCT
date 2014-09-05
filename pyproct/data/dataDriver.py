"""
Created on 6/8/2014

@author: victor
"""
from pyproct.data.matrix.matrixCalculator import MatrixCalculator
from pyproct.data.handler.dataHandler import DataHandler

class DataDriver(object):
    
    def __init__(self, parameters):
        data_handler = DataHandler(parameters)
        matrix_handler = MatrixCalculator.calculate(data_handler, 
                                                    parameters["matrix"])
        
        return data_handler, matrix_handler 
        
        
    
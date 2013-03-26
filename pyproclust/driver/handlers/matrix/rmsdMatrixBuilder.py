'''
Created on 13/02/2013

@author: victor
'''
from pyRMSD.RMSDCalculator import RMSDCalculator

class RMSDMatrixBuilder(object):

    def __init__(self):
        pass
    
    def build(self, trajectory_handler, matrix_creation_parameters):
        """
        Generates a matrix with the method used in the handler creation.
        
        @param trajectory_handler: 
        @param matrix_creation_parameters: 
        
        @return: The created matrix.
        """
        
        fit_selection = matrix_creation_parameters["fit_selection"]
        calc_selection = matrix_creation_parameters["calc_selection"]
        
        RMSDCalculator()
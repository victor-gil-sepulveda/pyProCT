"""
Created on 2/9/2014

@author: victor
"""
from pyRMSD.matrixHandler import MatrixHandler as pyRMSD_MatrixHandler

class LoaderMatrixCalculator(object):
    
    CALCULATION_METHOD = "matrix::load"
    
    def __init__(self, params):
        pass
    
    @classmethod
    def calculate(cls, data_handler, matrix_params):
        """
        :param matrix_params: The parameters to build the matrix. In this base case the only
        option is "load".
        
        Base parameters :
        
        {
            "method": STRING,
            "parameters":{
                ...
            }
        }
        
        Options:
        
        - "load": Load an already created matrix from disk
    
            "parameters":{
                "path": STRING
            }
        
        "path":  The path from where the matrix is going to be loaded.
        
        :return: A CondensedMatrix.
        """
        
        return pyRMSD_MatrixHandler.load_matrix(matrix_params["path"])

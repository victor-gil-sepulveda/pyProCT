"""
Created on 2/9/2014

@author: victor
"""
from pyRMSD.condensedMatrix import CondensedMatrix
import scipy.spatial.distance 

class EuclideanMatrixCalculator(object):
    
    CALCULATION_METHOD = "array::euclidean" 

    def __init__(self):
        pass
    
    @classmethod
    def calculate(cls, data_handler, matrix_params):
        return CondensedMatrix(scipy.spatial.distance.pdist(data_handler.get_data().get_all_elements(), 'euclidean'))
    
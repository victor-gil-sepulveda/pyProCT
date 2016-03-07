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
        data_array = data_handler.get_data().get_all_elements()
        if len(data_array.shape)==1: # 1 dimensional array
            distances = []
            for i in range(len(data_array)-1):
                for j in range(i+1, len(data_array)):
                    distances.append(abs(data_array[i]-data_array[j]))
            return CondensedMatrix(distances)
    
        else:
            return CondensedMatrix(scipy.spatial.distance.pdist(data_handler.get_data().get_all_elements(), 'euclidean'))
    
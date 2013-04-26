'''
Created on 07/02/2013

@author: victor
'''
import unittest
from pyproclust.driver.parameters import ProtocolParameters
from pyproclust.protocol.exploration.algorithmRunParametersGenerator import AlgorithmRunParametersGenerator

class MatrixHandlerMock:
    def __init__(self, matrix):
        self.distance_matrix = matrix

class MatrixMock:
    def __init__(self):
        self.row_length = 2000
    
    def calculateMean(self):
        return 2.5
    
    def calculateMax(self):
        return 4.0

class TestParameterGeneration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        parameters = ProtocolParameters.get_default_params("data/params.json")
        cls.parametersGenerator = AlgorithmRunParametersGenerator(parameters,MatrixHandlerMock(MatrixMock()))

    def test_get_gromos_parameters(self):
        expected = ([{'cutoff': 0.25}, 
                     {'cutoff': 0.5}, 
                     {'cutoff': 0.75}, 
                     {'cutoff': 1.0}, 
                     {'cutoff': 1.25}, 
                     {'cutoff': 1.5}, 
                     {'cutoff': 1.75}, 
                     {'cutoff': 2.0}, 
                     {'cutoff': 2.25}], [])
        
        self.assertItemsEqual(expected, TestParameterGeneration.parametersGenerator.get_gromos_parameters())
    
    def test_get_spectral_parameters(self):
        expected = ([{'k': 10, 'use_k_medoids': True}, 
                     {'k': 12, 'use_k_medoids': True}, 
                     {'k': 14, 'use_k_medoids': True}, 
                     {'k': 16, 'use_k_medoids': True}, 
                     {'k': 18, 'use_k_medoids': True}, 
                     {'k': 20, 'use_k_medoids': True}, 
                     {'k': 22, 'use_k_medoids': True}, 
                     {'k': 24, 'use_k_medoids': True}, 
                     {'k': 26, 'use_k_medoids': True}, 
                     {'k': 28, 'use_k_medoids': True},
                     {'k': 30, 'use_k_medoids': True}], [])

        self.assertItemsEqual(expected,TestParameterGeneration.parametersGenerator.get_spectral_parameters())
     
    def test_get_kmedoids_parameters(self):
        expected = ([{'seeding_type': 'GROMOS', 'k': 10, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 12, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 14, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 16, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 18, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 20, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 22, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 24, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 26, 'seeding_max_cutoff': 2.5}, 
                     {'seeding_type': 'GROMOS', 'k': 28, 'seeding_max_cutoff': 2.5},
                     {'seeding_type': 'GROMOS', 'k': 30, 'seeding_max_cutoff': 2.5}], [])

        self.assertItemsEqual(expected, TestParameterGeneration.parametersGenerator.get_kmedoids_parameters())
     
    def test_get_random_parameters(self):
        expected = ([{'num_clusters': 10}, 
                     {'num_clusters': 12}, 
                     {'num_clusters': 14}, 
                     {'num_clusters': 16}, 
                     {'num_clusters': 18}, 
                     {'num_clusters': 20}, 
                     {'num_clusters': 22}, 
                     {'num_clusters': 24}, 
                     {'num_clusters': 26}, 
                     {'num_clusters': 28},
                     {'num_clusters': 30}], [])

        self.assertItemsEqual(expected, TestParameterGeneration.parametersGenerator.get_random_parameters())
    
    ## get_hierarchical_parameters and get_dbscan_parameters depend on functions that may be tested apart.
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
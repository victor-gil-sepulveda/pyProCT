'''
Created on 07/02/2013

@author: victor
'''
import unittest
from pyproclust.driver.parameters import ProtocolParameters
from pyRMSD.condensedMatrix import CondensedMatrix
import pyproclust.algorithms.gromos.parametersGeneration as gromosParametersGeneration
import pyproclust.algorithms.kmedoids.parametersGeneration as kmedoidsParametersGeneration
import pyproclust.algorithms.random.parametersGeneration as randomParametersGeneration
import pyproclust.algorithms.spectral.parametersGeneration as spectralParametersGeneration

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
        cls.parameters = ProtocolParameters.get_default_params("data/params.json")
        distances = [94, 6, 43, 14, 96,
                        18, 59, 54, 69,
                            56, 96, 69,
                                54, 50,
                                     8]
        cls.matrix_1 = CondensedMatrix(distances)

    def test_get_gromos_parameters(self):
        expected = ([{'cutoff': 8.15},
                    {'cutoff': 16.3},
                    {'cutoff': 24.45},
                    {'cutoff': 32.6},
                    {'cutoff': 40.75},
                    {'cutoff': 48.9},
                    {'cutoff': 57.05},
                    {'cutoff': 65.2},
                    {'cutoff': 73.35},
                    {'cutoff': 81.5},
                    {'cutoff': 89.65}], [])
        
        parametersGenerator = gromosParametersGeneration.ParametersGenerator(self.parameters, 
                                                                             MatrixHandlerMock(self.matrix_1),
                                                                             10)
        parameters = parametersGenerator.get_parameters()[0]
        for i in  range(len(parameters)):
            self.assertAlmostEqual(parameters[i]["cutoff"], expected[0][i]["cutoff"]) 
    
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
        
        parametersGenerator = spectralParametersGeneration.ParametersGenerator(self.parameters, 
                                                                             MatrixHandlerMock(MatrixMock()))
        self.assertItemsEqual(expected,parametersGenerator.get_parameters())
     
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
        
        parametersGenerator = kmedoidsParametersGeneration.ParametersGenerator(self.parameters, 
                                                                             MatrixHandlerMock(MatrixMock()))
        self.assertItemsEqual(expected, parametersGenerator.get_parameters())
     
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

        parametersGenerator = randomParametersGeneration.ParametersGenerator(self.parameters, 
                                                                             MatrixHandlerMock(MatrixMock()))
        self.assertItemsEqual(expected, parametersGenerator.get_parameters())
    
    ## get_hierarchical_parameters and get_dbscan_parameters depend on functions that have been tested apart.
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
'''
Created on 27/05/2013

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.algorithms.gromos.parametersGeneration import ParametersGenerator

class MatrixHandlerMock:
    def __init__(self, matrix):
        self.distance_matrix = matrix

class Test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        distances = [94, 6, 43, 14, 96,
                        18, 59, 54, 69,
                            56, 96, 69,
                                54, 50,
                                     8]
        cls.matrix_1 = CondensedMatrix(distances)
        distances = [80, 1, 77, 73, 22, 
                        61, 48, 32, 31, 
                            73, 39, 19, 
                                90, 65, 
                                    89]
        cls.matrix_2 = CondensedMatrix(distances)
    
    def test_get_most_separated_elements(self):
        pg = ParametersGenerator({"clustering":{"algorithms":{"gromos":{"max":25}}}},MatrixHandlerMock(self.matrix_1))
        self.assertDictEqual(pg.get_most_separated_elements(), {'elements': (0, 5), 'value': 96.0})
        pg = ParametersGenerator({"clustering":{"algorithms":{"gromos":{"max":25}}}},MatrixHandlerMock(self.matrix_2))
        self.assertDictEqual(pg.get_most_separated_elements(), {'elements': (3, 4), 'value': 90.0})

    def test_get_most_separated_from_two_elements(self):
        pg = ParametersGenerator({"clustering":{"algorithms":{"gromos":{"max":25}}}},MatrixHandlerMock(self.matrix_1))
        self.assertDictEqual( pg.get_most_separated_from_two_elements(0, 5), {'mean_value': 81.5, 'element': 1, 'value': 94.0,})
        pg = ParametersGenerator({"clustering":{"algorithms":{"gromos":{"max":25}}}},MatrixHandlerMock(self.matrix_2))
        self.assertDictEqual( pg.get_most_separated_from_two_elements(3, 4),{'mean_value': 77.0, 'element': 5, 'value': 89.0})
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_get_most_separated_elementsName']
    unittest.main()
"""
Created on 13/11/2013

@author: victor
"""
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.tools.matrixTools import get_submatrix
import numpy

class Test(unittest.TestCase):


    def test_get_submatrix(self):
        old_matrix = CondensedMatrix([0.2, 1.,  0.3, 4.0, 6.1,
                                           0.5, 0.6, 0.7, 0.5,
                                                0.9, 0.8, 0.3,
                                                     0.4, 1.4,
                                                          2.9])
        expected = [0.9, 0.3,
                         1.4]

        numpy.testing.assert_array_almost_equal(expected, get_submatrix(old_matrix, [2,3,5]).get_data(),6)

        data = [1.0,  2.0,  3.0,   4.0,  5.0,  6.0,
                      7.0,  8.0,   9.0, 10.0, 11.0,
                            12.0, 13.0, 14.0, 15.0,
                                  16.0, 17.0, 18.0,
                                        19.0, 20.0,
                                              21.0]
        matrix = CondensedMatrix(data)

        all_elements_map = [2,4,6]
        new_matrix = get_submatrix(matrix, all_elements_map)
        expected1 = [13, 15,
                         20]
        numpy.testing.assert_array_equal(expected1, new_matrix.get_data())


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_get_submatrix']
    unittest.main()
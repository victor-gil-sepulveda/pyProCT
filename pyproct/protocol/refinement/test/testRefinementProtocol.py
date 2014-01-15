'''
Created on 17/09/2012

@author: victor
'''
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
import numpy
from pyproct.protocol.refinement.Refiner import Refiner
from pyproct.clustering.cluster import Cluster

class ClusterStub():
    def __init__(self,elements):
        self.all_elements = elements

class Test(unittest.TestCase):
    def test_get_map_and_submatrix(self):
        data = [1.0,  2.0,  3.0,   4.0,  5.0,  6.0,
                      7.0,  8.0,   9.0, 10.0, 11.0,
                            12.0, 13.0, 14.0, 15.0,
                                  16.0, 17.0, 18.0,
                                        19.0, 20.0,
                                              21.0]

        matrix = CondensedMatrix(data)
        submatrix, element_map = Refiner.get_cluster_submatrix_and_map(Cluster(None,[1,3,4]), matrix)
        self.assertItemsEqual([8., 9., 16.], submatrix.get_data())
        self.assertItemsEqual([1,3,4], element_map)

#     def test_redefine_cluster_with_map(self):
#

#     def test_recreate_matrix(self):
#         data = [1,  2,  3,  4,  5,  6,
#                     7,  8,  9, 10, 11,
#                        12, 13, 14, 15,
#                            16, 17, 18,
#                                19, 20,
#                                    21]
#         matrix = CondensedMatrix(data)
#
#         all_elements_map = [2,4,6]
#         new_matrix = recreate_matrix(matrix, len(all_elements_map), all_elements_map)
#         expected1 = [13, 15,
#                          20]
#         numpy.testing.assert_array_equal(expected1, new_matrix.get_data())
#
#         all_elements_map = [0,1,4,6]
#         new_matrix = recreate_matrix(matrix, len(all_elements_map), all_elements_map)
#         expected2 = [ 1, 4,  6,
#                          9, 11,
#                             20]
#         numpy.testing.assert_array_equal(expected2, new_matrix.get_data())


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_Name']
    unittest.main()
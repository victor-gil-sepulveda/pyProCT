'''
Created on 17/09/2012

@author: victor
'''
import unittest
from pyproct.protocol.refinementProtocol import redefine_clusters_with_map,\
    pick_all_elements_from_clusters, recreate_matrix
from pyRMSD.condensedMatrix import CondensedMatrix
import numpy

class ClusterStub():
    def __init__(self,elements):
        self.all_elements = elements

class Test(unittest.TestCase):

    def test_pick_all_elements_from_clusters(self):
        clusters = [ClusterStub([1,4,6,7,8]),ClusterStub([2,3,5,9,0])]
        expected_elements = [1,4,6,7,8,2,3,5,9,0]
        self.assertEqual(expected_elements, pick_all_elements_from_clusters(clusters)) 

    def test_redefine_clusters(self):
        expected_clusters = [ClusterStub([1,4,6,7,8]),ClusterStub([2,3,5,9,0])]
        clusters = [ClusterStub([0,1,2,3,4]),ClusterStub([5,6,7,8,9])]
        element_map = [1,4,6,7,8,2,3,5,9,0]
        new_clusters = redefine_clusters_with_map(clusters,element_map)
        trans_elements = pick_all_elements_from_clusters(new_clusters)
        expected_elems = pick_all_elements_from_clusters(expected_clusters)
        self.assertEqual(expected_elems, trans_elements) 
        
    def test_recreate_matrix(self):
        data = [1,  2,  3,  4,  5,  6,
                    7,  8,  9, 10, 11,
                       12, 13, 14, 15,
                           16, 17, 18,
                               19, 20,
                                   21]
        matrix = CondensedMatrix(data)
        
        all_elements_map = [2,4,6]
        new_matrix = recreate_matrix(matrix, len(all_elements_map), all_elements_map)
        expected1 = [13, 15,
                         20]
        numpy.testing.assert_array_equal(expected1, new_matrix.get_data())
        
        all_elements_map = [0,1,4,6]
        new_matrix = recreate_matrix(matrix, len(all_elements_map), all_elements_map)
        expected2 = [ 1, 4,  6, 
                         9, 11,
                            20]
        numpy.testing.assert_array_equal(expected2, new_matrix.get_data())
        
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_Name']
    unittest.main()
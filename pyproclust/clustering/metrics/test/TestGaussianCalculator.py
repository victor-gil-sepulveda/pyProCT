'''
Created on 13/06/2013

@author: victor
'''
import unittest

from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.clustering.metrics.test.data import  CH_table1
from pyproclust.clustering.cluster import Cluster
from pyproclust.clustering.clustering import Clustering
from pyproclust.clustering.metrics.gaussianSeparation import GaussianSeparationCalculator
from pyproclust.clustering.metrics.common import update_medoids
import numpy

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.matrix = CondensedMatrix(CH_table1)
         
        cls.clusterings = [Clustering([Cluster(None, [0,1,2,3]), Cluster(None, [4,5])]),
                            Clustering([Cluster(None, [0,1]), Cluster(None, [2,3]), Cluster(None, [4,5])])]
        
        update_medoids(cls.clusterings[0], cls.matrix)
        update_medoids(cls.clusterings[1], cls.matrix)

    def test_exponential_list_generation(self):
        numpy.testing.assert_almost_equal(GaussianSeparationCalculator.exponential_list_generation(self.clusterings[0], self.matrix),
                                          [2.61028153e-23], 5)
        numpy.testing.assert_almost_equal( GaussianSeparationCalculator.exponential_list_generation(self.clusterings[1], self.matrix),
                                           [7.78111377e-20, 4.78087732e-25, 2.93749177e-30], 5)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_']
    unittest.main()
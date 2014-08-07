"""
Created on 12/06/2013

@author: victor
"""
import unittest
from pyproct.clustering.metrics.test.data import squared_CH_table1
from pyproct.clustering.clustering import Clustering
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.metrics.DaviesBouldin import DaviesBouldinCalculator
from pyproct.clustering.metrics.common import update_medoids
from pyproct.clustering.cluster import Cluster
import numpy
class TestDaviesBouldin(unittest.TestCase):
   
    @classmethod
    def setUpClass(cls):
        cls.matrix = CondensedMatrix(squared_CH_table1)
         
        cls.clusterings = [Clustering([Cluster(None, [0,1,2,3]), Cluster(None, [4,5])]),
                            Clustering([Cluster(None, [0,1]), Cluster(None, [2,3]), Cluster(None, [4,5])])]
        update_medoids(cls.clusterings[0], cls.matrix)
        update_medoids(cls.clusterings[0], cls.matrix)
        
    def test_calculate_average_distances(self):
        
        self.assertItemsEqual(DaviesBouldinCalculator.calc_average_distances(self.clusterings[0], self.matrix),
                              [7,6])
       
    def test_max_db_term(self):
        numpy.testing.assert_almost_equal( DaviesBouldinCalculator.db_terms_for_cluster(0, [5.0, 6.0, 6.0], self.clusterings[1].clusters, self.matrix),
                             [1.0, 0.7857142857142857], 5)  
        numpy.testing.assert_almost_equal(  DaviesBouldinCalculator.db_terms_for_cluster(1, [5.0, 6.0, 6.0], self.clusterings[1].clusters, self.matrix),
                             [0.7058823529411765], 5)            
        numpy.testing.assert_almost_equal(  DaviesBouldinCalculator.db_terms_for_cluster(2, [5.0, 6.0, 6.0], self.clusterings[1].clusters, self.matrix),
                             [], 5)
        
        
    def test_db_eval(self):
        self.assertAlmostEqual(DaviesBouldinCalculator().evaluate(self.clusterings[1], self.matrix), 0.5686274509803922, 5)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
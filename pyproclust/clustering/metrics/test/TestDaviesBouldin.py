'''
Created on 12/06/2013

@author: victor
'''
import unittest
from pyproclust.clustering.metrics.test.data import squared_CH_table1
from pyproclust.clustering.clustering import Clustering
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.clustering.metrics.DaviesBouldin import DaviesBouldinCalculator
from pyproclust.clustering.metrics.common import update_medoids
from pyproclust.clustering.cluster import Cluster
import numpy
class TestDaviesBouldin(unittest.TestCase):
   
    @classmethod
    def setUpClass(cls):
        cls.matrix = CondensedMatrix(squared_CH_table1)
         
        cls.clusterings = [Clustering([Cluster(None, [0,1,2,3]), Cluster(None, [4,5])]),
                            Clustering([Cluster(None, [0,1]), Cluster(None, [2,3]), Cluster(None, [4,5])])]
 
     
    def test_calculate_average_distances(self):
        update_medoids(self.clusterings[0], self.matrix)
        self.assertItemsEqual(DaviesBouldinCalculator.calc_average_distances(self.clusterings[0], self.matrix),
                              [7,6])
       
    def test_max_db_term(self):
        update_medoids(self.clusterings[1], self.matrix)
        numpy.testing.assert_almost_equal( DaviesBouldinCalculator.db_terms_for_cluster(0, [5.0, 6.0, 6.0], self.clusterings[1].clusters, self.matrix),
                             [1.0, 0.7857142857142857], 5)  
        numpy.testing.assert_almost_equal(  DaviesBouldinCalculator.db_terms_for_cluster(1, [5.0, 6.0, 6.0], self.clusterings[1].clusters, self.matrix),
                             [0.7058823529411765], 5)            
        numpy.testing.assert_almost_equal(  DaviesBouldinCalculator.db_terms_for_cluster(2, [5.0, 6.0, 6.0], self.clusterings[1].clusters, self.matrix),
                             [], 5)
        
        
    def test_db_eval(self):
        update_medoids(self.clusterings[1], self.matrix)
        self.assertAlmostEqual(DaviesBouldinCalculator().evaluate(self.clusterings[1], self.matrix), 0.5686274509803922, 5)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
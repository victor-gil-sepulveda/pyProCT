'''
Created on 12/06/2013

@author: victor
'''
import unittest
from pyproclust.clustering.cluster import Cluster
from pyproclust.clustering.clustering import Clustering
from pyproclust.clustering.metrics.Dunn import DunnCalculator
from pyproclust.clustering.metrics.test.data import squared_CH_table1
from pyRMSD.condensedMatrix import CondensedMatrix


class Test(unittest.TestCase):
    
    def test_min_intracluster_distances(self):
        matrix = CondensedMatrix(squared_CH_table1)
        
        clusterings = [Clustering([Cluster(None, [0,1,2,3]), Cluster(None, [4,5])]),
                       Clustering([Cluster(None, [0,1]), Cluster(None, [2,3]), Cluster(None, [4,5])])] 
        
        expected = [5, 5]
        for i in range(len(clusterings)):
            self.assertEqual(DunnCalculator.min_intracluster_distances(clusterings[i], matrix), expected[i])
    
    def test_max_intercluster_distance(self):
        matrix = CondensedMatrix(squared_CH_table1)
        
        clusterings = [Clustering([Cluster(None, [0,1,2,3]), Cluster(None, [4,5])]),
                       Clustering([Cluster(None, [0,1]), Cluster(None, [2,3]), Cluster(None, [4,5])])] 
        
        expected = [21, 21]
        for i in range(len(clusterings)):
            self.assertEqual( DunnCalculator.max_intercluster_distance(clusterings[i], matrix), expected[i])
    
    def test_dunn_regression(self):
        matrix = CondensedMatrix(squared_CH_table1)
        
        clusterings = [Clustering([Cluster(None, [0,1,2,3]), Cluster(None, [4,5])]),
                       Clustering([Cluster(None, [0,1]), Cluster(None, [2,3]), Cluster(None, [4,5])])] 
        
        dunn_calculator = DunnCalculator()
        
        expected = [0.238095238095, 0.238095238095]
        
        for i in range(len(clusterings)):
            self.assertAlmostEqual( dunn_calculator.evaluate( clusterings[i], matrix), expected[i],6)
            
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
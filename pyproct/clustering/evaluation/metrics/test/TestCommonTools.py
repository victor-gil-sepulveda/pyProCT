"""
Created on 12/06/2013

@author: victor
"""
import unittest
from pyproct.clustering.evaluation.metrics.common import get_intra_cluster_distances,\
    get_inter_cluster_distances, get_distances_of_elements_to, update_medoids,\
    get_inter_cluster_prototype_distances
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.evaluation.metrics.test.data import CH_table1, squared_CH_table1
from pyproct.clustering.cluster import Cluster
import numpy
from pyproct.clustering.clustering import Clustering

class Test(unittest.TestCase):

    def test_get_intra_cluster_distances(self):
        matrix = CondensedMatrix(CH_table1)
        numpy.testing.assert_almost_equal(get_intra_cluster_distances(Cluster(None, [4,5]), matrix),[2.4494897427831779],5)
        numpy.testing.assert_almost_equal(get_intra_cluster_distances(Cluster(None, [1,3,5]), matrix),[2.4494897427831779, 3.8729833462074170, 3.8729833462074170],5)
        
        data = [0.0, 1.0, 1.0, 1.0,
                     1.0, 0.0, 0.0,
                          0.0, 0.0,
                               0.0]
        matrix = CondensedMatrix(data)
        expected_distance = 4
        self.assertEqual(expected_distance, numpy.sum(get_intra_cluster_distances(Cluster(None, range(5)), matrix)))

    def test_get_inter_cluster_distances(self):
        matrix = CondensedMatrix(squared_CH_table1)
        clusters = [Cluster(None, [1,2]),Cluster(None, [3,4]), Cluster(None, [5])]
        numpy.testing.assert_equal(get_inter_cluster_distances(0, 1, clusters, matrix), [6.0, 13.0, 6.0, 17.0])
        numpy.testing.assert_equal(get_inter_cluster_distances(1, 2, clusters, matrix), [15.0, 6.0])
        numpy.testing.assert_equal(get_inter_cluster_distances(0, 2, clusters, matrix), [15.0, 21.0])
    
    def test_get_inter_cluster_prototype_distances(self):
        matrix = CondensedMatrix(squared_CH_table1)
        clusters = [Cluster(2, [1,2]),Cluster(3, [3,4]), Cluster(0, [0,5])]
        numpy.testing.assert_equal( get_inter_cluster_prototype_distances(clusters, matrix), [6.0, 11.0, 11.0])
    
    def test_get_distances_of_elements_to(self):
        matrix = CondensedMatrix(list(squared_CH_table1))
        numpy.testing.assert_equal(get_distances_of_elements_to(3, [0,1,2,4,5], matrix), [11.0, 6.0, 6.0, 13.0, 15.0])
        
    def test_update_medois(self):
        clusters = [Cluster(None, [1,2]),Cluster(None, [3,4]), Cluster(None, [5])]
        clustering = Clustering(clusters)
        matrix = CondensedMatrix(squared_CH_table1)
        update_medoids(clustering, matrix)
        for c in clusters:
            self.assertNotEqual(c.prototype, None)
        
        self.assertItemsEqual([c.prototype for c in clusters], [1,3,5])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_get_intra_cluster_distances']
    unittest.main()
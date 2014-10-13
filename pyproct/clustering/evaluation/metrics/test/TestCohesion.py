"""
Created on 04/05/2012

@author: victor
"""
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.clustering import Clustering
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.evaluation.metrics.cython.cohesion import CohesionCalculator

class TestBoundedMetrics(unittest.TestCase):
    def test_cohesion(self):
        
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        
        clusters = [Cluster(None, elements=[0,1,2]),
                      Cluster(None, elements=[3,4])]
        clustering = Clustering(clusters)
        calculator = CohesionCalculator()
        
        # cohesion of cluster 1: 1/3 * 8 
        # cohesion of cluster 2: 1/2 * 10
        # max_cohesion =  11 (1/5 * 55)
        # final cohesion = 0.696945
        self.assertAlmostEqual(1-0.696945,calculator.evaluate(clustering, distances),places = 4)
        
    def test_cohesion_with_noise(self):
        # Element 2 is treated as noise
        distances =  CondensedMatrix([ 1., 2., 3., 4.,
                                           5., 6., 7.,
                                               8., 9.,
                                                  10.])
        
        clusters = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[3,4])]
        clustering = Clustering(clusters)
        calculator = CohesionCalculator()
        #[ 1.   3.   4.   
        #       6.   7.  
        #           10.]
        # cohesion of cluster 1: 1/2 
        # cohesion of cluster 2: 1/2 * 10
        # max_cohesion:  7.75 = (1/4 * 31)
        # final cohesion: 0.7096774193548387
        self.assertAlmostEqual(1-0.7096774193548387, calculator.evaluate(clustering, distances), places = 4)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_bounded_separation']
    unittest.main()
"""
Created on 02/04/2012

@author: victor
"""
import unittest
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.clustering import Clustering
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.evaluation.metrics.test import matrix
from pyproct.clustering.evaluation.metrics.meanMinimumDistance import MeanMinimumDistanceCalculator,\
    mean_function
from pyproct.clustering.algorithms.random import RandomAlgorithm
from pyproct.clustering.evaluation.metrics.cohesion import CohesionCalculator
from pyproct.clustering.evaluation.metrics.separation import SeparationCalculator

class TestMetrics(unittest.TestCase):
    
    def test_mean_function(self):
        self.assertEqual(mean_function([]),0)
        self.assertEqual(mean_function([1,2,3,4,5,6,7,8,9]), 5)
    
    def test_subsample(self):
        calculator = MeanMinimumDistanceCalculator(10)
        self.assertEqual( 14.5, calculator.subsample(10,20,  [4,6,3,24,7,12,9,17,20,43]))
        # If the percent gives less than 1 element, we use 1
        self.assertEqual( 43,  calculator.subsample(10,5,  [4,6,3,24,7,12,9,17,20,43]))
        # If the list has no elements, returns 0
        self.assertEqual( 0,  calculator.subsample(10,20,  []))
    
    def test_get_min_distances(self):
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters = [Cluster(None, elements=[0,1,2]),
                    Cluster(None, elements=[3,4])]
        
        calculator = MeanMinimumDistanceCalculator(10)
        
        min_dists, mean = calculator.get_mean_and_min_distances(clusters[0], clusters[1], distances)
        self.assertItemsEqual( [3.0, 6.0, 8.0], min_dists)
        self.assertAlmostEquals(12.33,mean,2)
        
    def test_get_distances_less_than_mean(self):
        calculator = MeanMinimumDistanceCalculator(10)
        self.assertItemsEqual(range(51), calculator.get_distances_less_than_mean(range(100),50))
    
    def test_subsampled_mean_min_dist(self):
        calculator = MeanMinimumDistanceCalculator(10)
        clusters = [Cluster(None, elements=[0,1,2]),
                    Cluster(None, elements=[3,4])]
        triangle = [ 1., 2., 3., 4., 
                         5., 6., 7., 
                             8., 9., 
                                10.]
        distances =  CondensedMatrix( triangle )
        self.assertEqual((8.0, 6.0),calculator.subsampled_mean_min_dist(clusters[0], clusters[1],  20, distances))
    
    def test_mini_evaluation(self):
        calculator = MeanMinimumDistanceCalculator(10)
        clusters = [Cluster(None, elements=[0,1,2]),
                    Cluster(None, elements=[3,4])]
        triangle = [ 1., 2., 3., 4., 
                         5., 6., 7., 
                             8., 9., 
                                10.]
        distances =  CondensedMatrix( triangle )
        clustering = Clustering(clusters)
        self.assertEqual(7.0, calculator.evaluate(clustering,distances,20))
        
    def test_full_run(self):
        condensed_matrix = CondensedMatrix(matrix)
        cmax = condensed_matrix.calculateMax()
        alg = RandomAlgorithm.RandomClusteringAlgorithm(condensed_matrix)
        values = []
        calculator =  MeanMinimumDistanceCalculator(10)
        for i in range(2,20):
            clustering = alg.perform_clustering({
                                                 "max_num_of_clusters":-1, 
                                                 "num_clusters":i
                                                 })
            values.append( calculator.evaluate(clustering, condensed_matrix, 30))
        self.assertTrue(max(values) < cmax)
        
    def test_cluster_cohesion_without_prototype(self):
        
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusters_2 = [Cluster(None, elements=[0,2,4]),
                      Cluster(None, elements=[1,3])]
        
        cohesion_calctor = CohesionCalculator()
        
        self.assertEqual(cohesion_calctor.evaluate_cluster(clusters_1[0],distances), 0.5)
        self.assertEqual(cohesion_calctor.evaluate_cluster(clusters_1[1],distances), 0.)
        self.assertEqual(cohesion_calctor.evaluate_cluster(clusters_1[2],distances), 5.0)
        self.assertEqual(cohesion_calctor.evaluate_cluster(clusters_2[0],distances), 5.0)
        self.assertEqual(cohesion_calctor.evaluate_cluster(clusters_2[1],distances), 3.0)

    def test_cluster_mixed_cohesion_wo_prot(self):
        
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusters_2 = [Cluster(None, elements=[0,2,4]),
                      Cluster(None, elements=[1,3])]
        
        sep_calctor = SeparationCalculator()
        
        self.assertEqual( sep_calctor._SeparationCalculator__between_cluster_distance(clusters_1[0],clusters_1[1],distances),7.0)
        self.assertEqual( sep_calctor._SeparationCalculator__between_cluster_distance(clusters_1[0],clusters_1[2],distances),20.0)
        self.assertEqual( sep_calctor._SeparationCalculator__between_cluster_distance(clusters_1[1],clusters_1[2],distances),17.0)
        self.assertEqual( sep_calctor._SeparationCalculator__between_cluster_distance(clusters_2[0],clusters_2[1],distances),34.0)
    
    def test_cluster_cohe_sep_wo_prot_eval(self):
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusters_2 = [Cluster(None, elements=[0,2,4]),
                      Cluster(None, elements=[1,3])]
        
        clusterization_1 = Clustering(clusters_1)
        clusterization_2 = Clustering(clusters_2)
        sep_calctor = SeparationCalculator()
        
        self.assertEqual(sep_calctor.cluster_separation(clusters_1[0],clusterization_1,1.,distances),27.0)
        self.assertEqual(sep_calctor.cluster_separation(clusters_1[1],clusterization_1,1.,distances),24.0)
        self.assertEqual(sep_calctor.cluster_separation(clusters_1[2],clusterization_1,1.,distances),37.0)
        self.assertEqual(sep_calctor.cluster_separation(clusters_2[0],clusterization_2,1.,distances),34.0)
        self.assertEqual(sep_calctor.cluster_separation(clusters_2[1],clusterization_2,1.,distances),34.0)
    
    def test_regression_separation_eval(self):
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters = [Cluster(None, elements=[0,1]),
                    Cluster(None, elements=[2]),
                    Cluster(None, elements=[3,4])]
        clustering = Clustering(clusters)
        
        sep_calctor = SeparationCalculator()
        self.assertEqual( sep_calctor.evaluate(clustering, distances,[1,1,1]), 27.0 + 24.0 + 37.0)
        self.assertEqual( sep_calctor.evaluate(clustering, distances), (1/0.5)*27.0 + (1/5.0)*37.0)
        
    def test_regression_cohesion_eval(self):
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters = [Cluster(None, elements=[0,1]),
                    Cluster(None, elements=[2]),
                    Cluster(None, elements=[3,4])]
        clustering = Clustering(clusters)
        
        cohesion_calctor = CohesionCalculator()
        self.assertEqual( cohesion_calctor.evaluate(clustering, distances), 5.5)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_cluster_cohesion']
    unittest.main()
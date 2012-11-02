'''
Created on 02/04/2012

@author: victor
'''
import unittest
from pyproclust.clustering.clusterization import Clustering
from pyproclust.clustering.cluster import Cluster
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
from pyproclust.clustering.metrics.clusteringMetrics import CohesionCalculator,\
    CohesionAndSeparationCalculator, mean_function,\
    MeanMinimumDistanceCalculator
from pyproclust.clustering.metrics.test import matrix
from pyproclust.algorithms.random import RandomAlgorithm

class Test(unittest.TestCase):
    
    def test_mean_function(self):
        self.assertEqual(mean_function([]),0)
        self.assertEqual(mean_function([1,2,3,4,5,6,7,8,9]), 5)
    
#    def test_mean_distance(self):
#        distances =  CondensedDistanceMatrix( [ 1., 2., 3., 4.,
#                                                    5., 6., 7., 
#                                                        8., 9., 
#                                                           10.])
#        clusters = [Cluster(None, elements=[0,1]),
#                    Cluster(None, elements=[2]),
#                    Cluster(None, elements=[3,4])]
#        
#        calculator = MeanMinimumDistanceCalculator()
#        
#        self.assertEqual(8.5, calculator.calculate_mean(clusters[1], clusters[2], distances))
#        self.assertEqual(5, calculator.calculate_mean(clusters[0], clusters[2], distances))
#        
    def test_subsample(self):
        calculator = MeanMinimumDistanceCalculator(10)
        self.assertEqual( 14.5, calculator.subsample(10,20,  [4,6,3,24,7,12,9,17,20,43]))
        # If the percent gives less than 1 element, we use 1
        self.assertEqual( 43,  calculator.subsample(10,5,  [4,6,3,24,7,12,9,17,20,43]))
        # If the list has no elements, returns 0
        self.assertEqual( 0,  calculator.subsample(10,20,  []))
    
    def test_get_min_distances(self):
        distances =  CondensedDistanceMatrix( [ 1., 2., 3., 4.,
                                                    5., 6., 7., 
                                                        8., 9., 
                                                           10.])
        clusters = [Cluster(None, elements=[0,1,2]),
                    Cluster(None, elements=[3,4])]
        
        calculator = MeanMinimumDistanceCalculator(10)
        
        min_dists, mean = calculator.get_mean_and_min_distances(clusters[0], clusters[1], distances)
        self.assertItemsEqual( [3.0, 6.0, 8.0], min_dists)
        self.assertAlmostEquals(6.166,mean,2)
        
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
        distances =  CondensedDistanceMatrix( triangle )
        self.assertItemsEqual((3.0, 6.0),calculator.subsampled_mean_min_dist(clusters[0], clusters[1],  20, distances))
    
    def test_mini_evaluation(self):
        calculator = MeanMinimumDistanceCalculator(10)
        clusters = [Cluster(None, elements=[0,1,2]),
                    Cluster(None, elements=[3,4])]
        triangle = [ 1., 2., 3., 4., 
                         5., 6., 7., 
                             8., 9., 
                                10.]
        distances =  CondensedDistanceMatrix( triangle )
        clustering = Clustering(clusters)
        self.assertEqual(4.5, calculator.evaluate(clustering,20,distances))
        
    def test_full_run(self):
        condensed_matrix = CondensedDistanceMatrix(matrix)
        cmin,cmax = condensed_matrix.get_minimum_and_maximum()
        del cmin
        alg = RandomAlgorithm.RandomClusteringAlgorithm(condensed_matrix)
        values = []
        calculator =  MeanMinimumDistanceCalculator(10)
        for i in range(0,20):
            clustering = alg.perform_clustering({"max_num_of_clusters":0, "num_clusters":i})
            values.append( calculator.evaluate(clustering, 30, condensed_matrix))
        self.assertTrue(max(values) < cmax)
        
    def test_cluster_cohesion_without_prototype(self):
        
        distances =  CondensedDistanceMatrix( [ 1., 2., 3., 4.,
                                                    5., 6., 7., 
                                                        8., 9., 
                                                           10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusters_2 = [Cluster(None, elements=[0,2,4]),
                      Cluster(None, elements=[1,3])]
        
        cohesion_calctor = CohesionCalculator()
        
        self.assertEqual(cohesion_calctor.evaluate(clusters_1[0],distances), 0.5)
        self.assertEqual(cohesion_calctor.evaluate(clusters_1[1],distances), 0.)
        self.assertEqual(cohesion_calctor.evaluate(clusters_1[2],distances), 5.0)
        self.assertEqual(cohesion_calctor.evaluate(clusters_2[0],distances), 5.0)
        self.assertEqual(cohesion_calctor.evaluate(clusters_2[1],distances), 3.0)

    def test_cluster_mixed_cohesion_wo_prot(self):
        
        distances =  CondensedDistanceMatrix( [ 1., 2., 3., 4.,
                                                    5., 6., 7., 
                                                        8., 9., 
                                                           10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusters_2 = [Cluster(None, elements=[0,2,4]),
                      Cluster(None, elements=[1,3])]
        
        cosep_calctor = CohesionAndSeparationCalculator()
        
        self.assertEqual( cosep_calctor._CohesionAndSeparationCalculator__clusters_mixed_cohesion_wo_prot(clusters_1[0],clusters_1[1],distances),7.0)
        self.assertEqual( cosep_calctor._CohesionAndSeparationCalculator__clusters_mixed_cohesion_wo_prot(clusters_1[0],clusters_1[2],distances),20.0)
        self.assertEqual( cosep_calctor._CohesionAndSeparationCalculator__clusters_mixed_cohesion_wo_prot(clusters_1[1],clusters_1[2],distances),17.0)
        self.assertEqual( cosep_calctor._CohesionAndSeparationCalculator__clusters_mixed_cohesion_wo_prot(clusters_2[0],clusters_2[1],distances),34.0)
    
    def test_cluster_cohe_sep_wo_prot_eval(self):
        distances =  CondensedDistanceMatrix( [ 1., 2., 3., 4.,
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
        cosep_calctor = CohesionAndSeparationCalculator()
        
        self.assertEqual(cosep_calctor._CohesionAndSeparationCalculator__noproto_eval(clusters_1[0],clusterization_1,1.,distances),27.0)
        self.assertEqual(cosep_calctor._CohesionAndSeparationCalculator__noproto_eval(clusters_1[1],clusterization_1,1.,distances),24.0)
        self.assertEqual(cosep_calctor._CohesionAndSeparationCalculator__noproto_eval(clusters_1[2],clusterization_1,1.,distances),37.0)
        self.assertEqual(cosep_calctor._CohesionAndSeparationCalculator__noproto_eval(clusters_2[0],clusterization_2,1.,distances),34.0)
        self.assertEqual(cosep_calctor._CohesionAndSeparationCalculator__noproto_eval(clusters_2[1],clusterization_2,1.,distances),34.0)
    
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_cluster_cohesion']
    unittest.main()
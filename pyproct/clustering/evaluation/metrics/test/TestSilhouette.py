"""
Created on 27/01/2014

@author: victor
"""
import unittest

from pyproct.clustering.clustering import Clustering
from pyproct.clustering.cluster import Cluster

import pyproct.clustering.evaluation.metrics.test.data.example_clustering_1 as data
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.evaluation.metrics.cython.silhouette import SilhouetteCoefficientCalculator
from pyproct.data.matrix.matrixHandler import MatrixHandler


class Test(unittest.TestCase):

    def testSilhouetteSpecialCase(self):
        clustering = Clustering.from_dic(data.clustering_01)
        mh = MatrixHandler({
                                "method": "load",
                                "parameters":{
                                    "path": "data/example_clustering_1_matrix"
                                }
                            }
        )
        s = SilhouetteCoefficientCalculator()
        matrix =  mh.create_matrix(None)
        print s.evaluate(clustering, matrix)

    def test_get_average_distance(self):
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusters_2 = [Cluster(None, elements=[0,2,4]),
                      Cluster(None, elements=[1,3])]
        
        sil_calc = SilhouetteCoefficientCalculator()
        self.assertEqual( sil_calc._SilhouetteCoefficientCalculator__get_average_distance_with_cluster(0,clusters_1[0],distances),0.5)
        self.assertEqual( sil_calc._SilhouetteCoefficientCalculator__get_average_distance_with_cluster(1,clusters_1[0],distances),0.5)
        self.assertEqual( sil_calc._SilhouetteCoefficientCalculator__get_average_distance_with_cluster(2,clusters_1[1],distances),0.0)
        self.assertEqual( sil_calc._SilhouetteCoefficientCalculator__get_average_distance_with_cluster(2,clusters_1[2],distances),8.5)

        self.assertEqual( sil_calc._SilhouetteCoefficientCalculator__get_average_distance_with_cluster(0,clusters_2[0],distances),2.0)
        self.assertEqual( sil_calc._SilhouetteCoefficientCalculator__get_average_distance_with_cluster(0,clusters_2[1],distances),2.0)
        
    
    def test_one_element_silhouette(self):
        
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
        
        sil_calc = SilhouetteCoefficientCalculator()
        
        self.assertEqual(  sil_calc._SilhouetteCoefficientCalculator__one_element_silhouette(0,clusters_1[0],clusterization_1,distances),0.5)
        self.assertEqual(  sil_calc._SilhouetteCoefficientCalculator__one_element_silhouette(1,clusters_1[0],clusterization_1,distances),0.8)
        self.assertAlmostEqual(sil_calc._SilhouetteCoefficientCalculator__one_element_silhouette(2,clusters_1[1],clusterization_1,distances),0.714, places = 3)
        self.assertAlmostEqual( sil_calc._SilhouetteCoefficientCalculator__one_element_silhouette(3,clusters_1[2],clusterization_1,distances),-0.55)
        self.assertAlmostEqual( sil_calc._SilhouetteCoefficientCalculator__one_element_silhouette(0,clusters_2[0],clusterization_2,distances),-0.333, places = 3)
        self.assertAlmostEqual( sil_calc._SilhouetteCoefficientCalculator__one_element_silhouette(1,clusters_2[1],clusterization_2,distances),-0.2777, places = 3)

    def test_one_cluster_silhouette(self):
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusterization_1 = Clustering(clusters_1)
        
        sil_calc = SilhouetteCoefficientCalculator()
        
        self.assertItemsEqual( sil_calc._SilhouetteCoefficientCalculator__one_cluster_partial_silhouette(clusters_1[0],clusterization_1,distances),[0.5, 0.80000000000000004])
        self.assertItemsEqual( sil_calc._SilhouetteCoefficientCalculator__one_cluster_partial_silhouette(clusters_1[1],clusterization_1,distances),[0.7142857142857143])
        self.assertItemsEqual( sil_calc._SilhouetteCoefficientCalculator__one_cluster_partial_silhouette(clusters_1[2],clusterization_1,distances),[-0.55000000000000004, -0.45000000000000001])

    def test_one_clusterization_silhouette(self):
        distances =  CondensedMatrix( [ 1., 2., 3., 4.,
                                            5., 6., 7., 
                                                8., 9., 
                                                   10.])
        clusters_1 = [Cluster(None, elements=[0,1]),
                      Cluster(None, elements=[2] ),
                      Cluster(None, elements=[3,4])]
        
        clusterization_1 = Clustering(clusters_1)
        sil_calc = SilhouetteCoefficientCalculator()
        expected = [0.5, 0.80000000000000004, -0.55000000000000004, -0.45000000000000001, 0.7142857142857143]
        
        self.assertItemsEqual(sil_calc._SilhouetteCoefficientCalculator__one_clusterization_partial_silhouette(clusterization_1,distances),expected)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testSilhouetteSpecialCase']
    unittest.main()
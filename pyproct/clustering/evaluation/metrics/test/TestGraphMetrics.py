"""
Created on 13/08/2012

@author: victor
"""
import unittest
from pyproct.clustering.clustering import Clustering
from pyproct.clustering.cluster import Cluster
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.evaluation.metrics.cython.graph.tools import get_cluster_and_complementary

class TestGraphMetrics(unittest.TestCase):

    def test_getClusterAndComplementary(self):
        clustering = Clustering([Cluster(1,range(5)),Cluster(5,range(5,10)),Cluster(10,range(10,20))])
        A,Acomp = get_cluster_and_complementary(1, clustering.clusters)
        A.sort() 
        Acomp.sort()
        self.assertItemsEqual(A, [0, 1, 2, 3, 4] )
        self.assertItemsEqual(Acomp,[5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    
#     def test_W(self):
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         self.assertEqual(W(range(3),range(3,6),matrix), 37)
#         
#     def test_d(self):
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         self.assertEqual(d(2, matrix),13.0)
#         
#     def test_vol(self):
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         
#         self.assertEqual( vol([2,4,3],matrix),71)
#         
#     def test_all_cut(self):
#         clustering = Clustering([Cluster(1,[0,1]),Cluster(2,[2,3]),Cluster(4,[4,5])])
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         
#         self.assertEqual(  all_cut(clustering,matrix), 52)
#         
#     def test_single_cut(self):
#         clustering = Clustering([Cluster(1,[0,1]),Cluster(2,[2,3]),Cluster(4,[4,5])])
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         
#         self.assertEqual(  all_cut(clustering,matrix), 52)
#         
#     def test_regression_MinMaxCut(self):
#         clustering = Clustering([Cluster(1,[0,1]),Cluster(2,[2,3]),Cluster(4,[4,5])])
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         
#         calculator = MinMaxCut()
#         
#         self.assertAlmostEqual(calculator.evaluate(clustering, matrix),6.34375,4)
#     
#     def test_regression_NCut(self):
#         clustering = Clustering([Cluster(1,[0,1]),Cluster(2,[2,3]),Cluster(4,[4,5])])
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         
#         calculator = NCut()
#         
#         self.assertAlmostEqual( calculator.evaluate(clustering, matrix), 7.4983003265331849, 4)
#         
#     def test_regression_RatioCut(self):
#         clustering = Clustering([Cluster(1,[0,1]),Cluster(2,[2,3]),Cluster(4,[4,5])])
#         matrix_data = [1., 4., 6., 2., 5.,
#                            3., 9., 7., 2.,
#                                4., 1., 1.,
#                                    9., 3.,
#                                        8.]
#         matrix = CondensedMatrix(matrix_data)
#         
#         calculator = MinMaxCut()
#         
#         self.assertAlmostEqual( calculator.evaluate(clustering, matrix),6.34375,4)
#         
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_getClusterAndComplementary']
    unittest.main()
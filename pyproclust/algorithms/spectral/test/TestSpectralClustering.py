'''
Created on 09/01/2013

@author: victor
'''
import unittest
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm
from scipy.spatial.distance import pdist

class TestSpectralClustering(unittest.TestCase):
    
    def test_calculate_adjacency_matrix_regression(self):
        
        data = [1., 4., 6., 2., 5.,
                    3., 9., 7., 2.,
                        4., 1., 1.,
                             9.,3.,
                                8.]
        
        expected_W = [[  1.0, 1.35375977e-01, 0.0, 0.0,  3.35454941e-04, 0.0],
                     [  1.35375977e-01, 1.0, 0.0, 0.0, 0.0, 3.35454941e-04],
                     [  0.0, 0.0, 1.0, 0.0, 1.35375977e-01, 1.35375977e-01],
                     [  0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                     [  3.35454941e-04, 0.0, 1.35375977e-01, 0.0, 1.0, 0.0],
                     [  0.0, 3.35454941e-04, 1.35375977e-01, 0.0, 0.0, 1.0]]
        matrix  = CondensedMatrix(data)
   
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix, sigma_sq = 0.5)
        numpy.testing.assert_almost_equal(W,expected_W)
           
    def test_calculate_laplacian(self):
        data = [1., 4., 6., 2., 5.,
                    3., 9., 7., 2.,
                        4., 1., 1.,
                             9.,3.,
                                8.]
          
        expected_regression = [[ 0.11949462, -0.11919928, 0., 0.,  -0.00029537, 0. ],
                                 [-0.11919928, 0.11949462, 0., 0., 0.,  -0.00029537],
                                 [ 0., 0., 0.21306437, 0.,  -0.10653218, -0.10653218],
                                 [ 0., 0., 0., 0., 0., 0. ],
                                 [-0.00029537, 0., -0.11919928, 0., 0.11949462, 0.],
                                 [ 0., -0.00029537, -0.11919928, 0., 0., 0.11949462]]
        
        matrix  = CondensedMatrix(data)
   
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix, sigma_sq = 0.5)
           
        L_PYTHON, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "PYTHON")
           
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix, sigma_sq = 0.5)
           
        L_NUMPY, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NUMPY")
           
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix, sigma_sq = 0.5)
           
        L_NUMPY_PURE, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NUMPY_PURE")
        
        numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY, 3)
        numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY_PURE, 3)
        numpy.testing.assert_almost_equal(L_PYTHON, expected_regression, 3)
          
    def test_calculate_degree_matrix(self):
        W = [[1., 4., 6.],
            [4., 1., 7.],
            [6., 7., 1.]]
         
        expected_D = [11., 12., 14.]
        D = SpectralClusteringAlgorithm.calculate_degree_matrix(W)
        numpy.testing.assert_almost_equal(D, expected_D, 8)
        
        expected_D_inv = [1/11., 1/12., 1/14.]
        D = SpectralClusteringAlgorithm.calculate_inverse_degree_matrix(W)
        numpy.testing.assert_almost_equal(D, expected_D_inv, 8)
         
    def test_naive_case_1(self):
#         1       5         8
#         |       |         |
#         0 - 3   4     6 - 7 
#         |                 |
#         2                 9
        points = [(0,0),(0,1),(0,-1),(1,0),
                  (3,0),(3,1),
                  (6,0),(7,0),(7,1),(7,-1)]
        
        matrix = CondensedMatrix(pdist(points))
        s_algo = SpectralClusteringAlgorithm(matrix, max_clusters = 3,  verbose = True)
        clusters = s_algo.perform_clustering({"k":3}).clusters
        
        for c in clusters:
            print c
            self.assertIn(c.all_elements, [[0, 1, 2, 3],[6, 7, 8, 9],[4, 5]])

    def test_naive_case_2(self):
#                     2 
#                     |
#         0 - 4       3
#         |       
#         1       
        points = [(0,0), (0,-1), (6,1), (6,0), (1,0)]
        
        matrix = CondensedMatrix(pdist(points))
        s_algo = SpectralClusteringAlgorithm(matrix, max_clusters = 2, sigma_sq = 0.5, verbose = True)
        clusters = s_algo.perform_clustering({"k":2}).clusters
        
        for c in clusters:
            print c
            self.assertIn(c.all_elements, [[0,1,4],[2,3]])
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
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
         
        expected_W = [[ 1.0, 2.14080811e-02, 0.0, 0.0, 2.38418579e-07, 0.0], 
                      [ 2.14080811e-02, 1.0, 0.0, 0.0, 0.0, 2.38418579e-07], 
                      [ 0.0, 0.0, 1.0, 0.0, 2.14080811e-02, 2.14080811e-02], 
                      [ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0], 
                      [ 2.38418579e-07, 0.0, 2.14080811e-02, 0.0, 1.0, 0.0], 
                      [ 0.0, 2.38418579e-07, 2.14080811e-02, 0.0, 0.0, 1.0]]
         
        matrix  = CondensedMatrix(data)
  
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix)
        numpy.testing.assert_almost_equal(W,expected_W)
          
    def test_calculate_laplacian(self):
        data = [1., 4., 6., 2., 5.,
                    3., 9., 7., 2.,
                        4., 1., 1.,
                             9.,3.,
                                8.]
         
        expected_regression = [[0.9444444440305233, -0.0011893378105014563, 0.0, 0.0, -1.3245476715439963e-08, 0.0],
                               [-0.0009730946039780974, 0.9545454531908035, 0.0, 0.0, 0.0, -1.0837208463954084e-08],
                               [0.0, 0.0, 0.9230769202113152, 0.0, -0.0016467755194753408, -0.0016467755194753408],
                               [0.0, 0.0, 0.0, 0.9677419364452362, 0.0, 0.0],
                               [-8.830317810293309e-09, 0.0, -0.0007928918930701911, 0.0, 0.9629629626870155, 0.0],
                               [0.0, -1.2548346361995755e-08, -0.0011267410591244698, 0.0, 0.0, 0.9473684206604958]]
 
        matrix  = CondensedMatrix(data)
  
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix)
          
        L_PYTHON, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "PYTHON")
          
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix)
          
        L_NUMPY, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NUMPY")
          
        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix)
          
        L_NUMPY_PURE, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NUMPY_PURE")
  
        numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY, 3)
        numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY_PURE, 3)
        numpy.testing.assert_almost_equal(L_PYTHON, expected_regression, 3)
         
          
    def test_calculate_degree_matrix(self):
        data = [1., 4., 6., 2., 5.,
                    3., 9., 7., 2.,
                        4., 1., 1.,
                             9.,3.,
                                8.]
        expected_D = [18., 22., 13., 31., 27., 19.]
        matrix  = CondensedMatrix(list(data))
        D = SpectralClusteringAlgorithm.calculate_degree_matrix(matrix)
        numpy.testing.assert_almost_equal(D, expected_D, 8)
        
    def test_naive_case(self):
#         1       5         8
#         |       |         |
#         0 - 3   4     6 - 7 
#         |                 |
#         2                 9
        points = [(0,0),(0,1),(0,-1),(1,0),
                  (3,0),(3,1),
                  (6,0),(7,0),(7,1),(7,-1)]
        matrix = CondensedMatrix(pdist(points))
        s_algo = SpectralClusteringAlgorithm(matrix, 3)
        clusters = s_algo.perform_clustering({"k":3}).clusters
        
        for c in clusters:
            self.assertIn(c.all_elements, [[0, 1, 2, 3],[6, 7, 8, 9],[4, 5]])
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
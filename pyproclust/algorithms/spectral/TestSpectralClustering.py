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
        expected_W = [[1.,0.87792969,0.59423828,0.45825195,0.77099609,0.52197266],
                    [0.87792969,1.,0.67675781,0.31005859,0.40234375,0.77099609],
                    [0.59423828,0.67675781,1.,0.59423828,0.87792969,0.87792969],
                    [0.45825195,0.31005859,0.59423828,1.,0.31005859,0.67675781],
                    [0.77099609,0.40234375,0.87792969,0.31005859,1.,0.35302734],
                    [0.52197266,0.77099609,0.87792969,0.67675781,0.35302734,1.]]
        matrix  = CondensedMatrix(data)

        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix)
        numpy.testing.assert_almost_equal(W,expected_W)
        
    def test_calculate_adjacency_matrix(self):
        data = [1., 4., 6., 2., 5.,
                    3., 9., 7., 2.,
                        4., 1., 1.,
                             9.,3.,
                                8.]
        expected_regression = [[0.94433594,-0.04876709,-0.03302002,-0.02545166,-0.04284668,-0.0289917],
                             [-0.03991699, 0.95458984,-0.03076172,-0.01409149,-0.01829529,-0.03503418],
                             [-0.04571533,-0.05206299,0.92285156,-0.04571533,-0.06750488,-0.06750488],
                             [-0.01478577,-0.01000214,-0.01916504,0.96777344,-0.01000214,-0.02183533],
                             [-0.02854919,-0.01490021,-0.03250122,-0.01148224,0.96289062,-0.01307678],
                             [-0.02746582,-0.04058838,-0.04620361,-0.03561401,-0.01858521,0.94726562]]
        
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
            print c
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
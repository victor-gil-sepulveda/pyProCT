'''
Created on 09/01/2013

@author: victor
'''
import unittest
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm
from scipy.spatial.distance import pdist

class TestSpectralClustering(unittest.TestCase):



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

        L_PYTHON, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NORM_PYTHON")

#         W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix, sigma_sq = 0.5)
#
#         L_NUMPY, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NORM_NUMPY")

        W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix, sigma_sq = 0.5)

        L_NUMPY_PURE, D = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NORM_NUMPY_PURE")

#         numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY, 3)
        numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY_PURE, 3)
        numpy.testing.assert_almost_equal(L_PYTHON, expected_regression, 3)


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
            self.assertIn(c.all_elements, [[0,1,4],[2,3]])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
"""
Created on 09/01/2013

@author: victor
"""
import unittest
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm
from scipy.spatial.distance import pdist

class TestSpectralClustering(unittest.TestCase):

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
        s_algo = SpectralClusteringAlgorithm(matrix, sigma = 20, max_clusters = 3, use_k_medoids = False, verbose = True)
        clusters = s_algo.perform_clustering({"k":3}).clusters

        # sometimes works, sometimes not (due to kmeans/medoids unstabilities)
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
        s_algo = SpectralClusteringAlgorithm(matrix, max_clusters = 2, sigma_sq = 0.5, use_k_medoids = False, verbose = True)
        clusters = s_algo.perform_clustering({"k":2}).clusters

        for c in clusters:
            self.assertIn(c.all_elements, [[0,1,4],[2,3]])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()

"""
Created on 07/05/2014

@author: victor
"""
import unittest
import numpy

import pyproct.clustering.algorithms.spectral.cython.spectralTools as SpectralTools
from pyRMSD.condensedMatrix import CondensedMatrix

class Test(unittest.TestCase):

    def test_order_by_eigenvalue(self):
        vectors = numpy.array([[ 1, 2, 3],
                               [ 4, 5, 6],
                               [ 7, 8,9 ]])

        values = numpy.array([3,1,2])

        expected_vectors = numpy.array( [[2, 3, 1,],
                                         [5, 6, 4,],
                                         [8, 9, 7,]])

        v = SpectralTools.order_by_eigenvalues(vectors, values)

        numpy.testing.assert_array_equal(expected_vectors, v)

        values = numpy.array([3,1])

        v = SpectralTools.order_by_eigenvalues(vectors, values)

        expected_vectors = numpy.array( [[2, 3],
                                         [5, 6],
                                         [8, 9]])

        v = SpectralTools.order_by_eigenvalues(vectors, values)

    def test_calculate_degree_matrix(self):
        W_data = numpy.array([4., 6.,
                                  7.])
        expected_D = [10., 11., 13.]
        D = SpectralTools.calculate_degree_matrix(CondensedMatrix(W_data))
        numpy.testing.assert_almost_equal(D, expected_D, 8)

    def test_calculate_adjacency_matrix(self):
        data = numpy.array([1., 4., 6., 2., 5.,
                                3., 9., 7., 2.,
                                    4., 1., 1.,
                                         9.,3.,
                                            8.])
        matrix  = CondensedMatrix(data)

        expected_W = numpy.array(    [ 0.13537598,  0.        ,  0.        ,  0.00033545,  0.        ,
                                                    0.        ,  0.        ,  0.        ,  0.00033545,
                                                                 0.        ,  0.13537598,  0.13537598,
                                                                              0.        ,  0.        ,
                                                                                           0.        ])

        W = SpectralTools.calculate_fully_connected_adjacency_matrix(matrix, sigma_sq = 0.5)
        numpy.testing.assert_almost_equal(W.get_data(),expected_W,3)# float comparison

    def test_calculate_adjacency_matrix_with_sigma_estimation_REGRESSION(self):
        sigmas = [1.,2.,3.,4.,5.,6.]

        data = numpy.array([1., 4., 6., 2., 5.,
                                3., 9., 7., 2.,
                                    4., 1., 1.,
                                         9.,3.,
                                            8.])
        matrix  = CondensedMatrix(data)

        expected_W = numpy.array([ 6.06530660e-01,  4.82794999e-03,  1.23409804e-04,  4.49328964e-01,  1.55038536e-02,
                                                    2.23130160e-01,  4.00652974e-05,  7.44658307e-03,  7.16531311e-01,
                                                                     2.63597138e-01,  9.35506985e-01,  9.45959469e-01,
                                                                                      1.74223746e-02,  6.87289279e-01,
                                                                                                       1.18441829e-01])

        numpy.testing.assert_almost_equal(expected_W,\
                                          SpectralTools.calculate_fully_connected_adjacency_matrix_with_sigma_estimation(matrix, sigmas).get_data(),\
                                          8)

    def test_force_sparsity(self):
        data = numpy.array([1., 4., 6., 2., 5.,
                                3., 9., 7., 2.,
                                    4., 1., 1.,
                                         9.,3.,
                                            8.])
        W  = CondensedMatrix(data)
        SpectralTools.force_sparsity(W)
        W_expected = [ 0., 0., 6., 0., 5., 0., 9., 7., 0., 0., 0., 0., 9., 0., 8.]
        numpy.testing.assert_array_equal(W_expected, W.get_data())


    def test_calculate_unnormalized_laplacian(self):
        data = numpy.array([1., 4., 6.,
                                3., 9.,
                                    4.])
        D = numpy.array([2,2,2,2])
        W  = CondensedMatrix(data)

        L_expected = [[1., -1., -4., -6.],
                      [-1,  1., -3., -9.],
                      [-4, -3.,  1., -4.],
                      [-6, -9., -4., 1.]]

        L = SpectralTools.calculateUnnormalizedLaplacian(W,D)
        numpy.testing.assert_array_equal(L_expected, L )

    def test_eigencalculation(self):
        L = numpy.array([  [1., 0., 0.,  0.],
                          [ 0,  1.,  0., -4.],
                          [ 0,  0.,  1., 0.],
                          [ 0, -4.,  0.,  1.]])
        max_clusters = 2
        v  = SpectralTools.calculateUnnormalizedEigenvectors(L, max_clusters, False)

        L = numpy.array([  [1., 0., 0.,  0.],
                          [ 0,  1.,  0., -4.],
                          [ 0,  0.,  1., 0.],
                          [ 0, -4.,  0.,  1.]])
        v  = SpectralTools.calculateUnnormalizedEigenvectors(L, max_clusters, True)

        D = [1.,2.,3.,4.]

        L = numpy.array([  [1., 0., 0.,  0.],
                          [ 0,  1.,  0., -4.],
                          [ 0,  0.,  1., 0.],
                          [ 0, -4.,  0.,  1.]])
        v  = SpectralTools.calculateNormalizedEigenvectors(L, D, max_clusters, False)

        L = numpy.array([  [1., 0., 0.,  0.],
                          [ 0,  1.,  0., -4.],
                          [ 0,  0.,  1., 0.],
                          [ 0, -4.,  0.,  1.]])
        v  = SpectralTools.calculateNormalizedEigenvectors(L, D, max_clusters, True)
        self.fail("CHECK EIGENCALCULATIONS AGAIN, VALUES FOR THIS METHODS ARE DIFFERENT")
        
    def test_calculate_degree_matrix_2(self):
        """
        Test provided by Nancy-Sarah Yacovzada.
        """
        
#         0--3
#         |\  \ 
#         5 4  2
#         
#         X = np.array(
#              [[0.0, 0.0, 1.0, 1.0, 1.0],
#               [0.0, 0.0, 1.0, 0.0, 0.0],
#               [1.0, 1.0, 0.0, 0.0, 0.0],
#               [1.0, 0.0, 0.0, 0.0, 0.0],
#               [1.0, 0.0, 0.0, 0.0, 0.0]]
#               )
        data = [0.0, 1.0, 1.0, 1.0,
                     1.0, 0.0, 0.0,
                          0.0, 0.0,
                               0.0]
        
        W = CondensedMatrix(data)
        
        self.assertListEqual([ 3., 1., 2., 1., 1.], SpectralTools.calculate_degree_matrix(W))
        
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_order_by_eigenvalue']
    unittest.main()
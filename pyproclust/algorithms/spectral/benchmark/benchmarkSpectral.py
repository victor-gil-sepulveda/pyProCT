'''
Created on 14/01/2013

@author: victor
'''
import numpy
import time 
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm

if __name__ == '__main__':
    MAX_ELEMENTS = 200
    DATA_LEN = (MAX_ELEMENTS * (MAX_ELEMENTS-1))/2
    numpy_matrix_data = numpy.abs(numpy.random.rand(DATA_LEN))
    matrix  = CondensedMatrix(numpy_matrix_data)

    W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix)

    t1 = time.time()
    L_PYTHON = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "PYTHON")
    t2 = time.time()
    print 'Calculating took %0.3f s' % (t2-t1)
     
    t1 = time.time()
    L_NUMPY = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NUMPY")
    t2 = time.time()
    print 'Calculating took %0.3f s' % (t2-t1)
     
    t1 = time.time()
    L_NUMPY_PURE = SpectralClusteringAlgorithm.calculate_laplacian(W, matrix, "NUMPY_PURE")
    t2 = time.time()
    print 'Calculating took %0.3f s' % (t2-t1)
    
    numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY, 3)
    numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY_PURE, 3)
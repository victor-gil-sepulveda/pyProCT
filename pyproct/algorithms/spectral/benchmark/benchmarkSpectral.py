"""
Created on 14/01/2013

@author: victor
"""
import numpy
import time 
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm

if __name__ == '__main__':
    MAX_ELEMENTS = 2000
    DATA_LEN = (MAX_ELEMENTS * (MAX_ELEMENTS-1))/2
    numpy_matrix_data = numpy.abs(numpy.random.rand(DATA_LEN))
    matrix  = CondensedMatrix(numpy_matrix_data)

    W = SpectralClusteringAlgorithm.calculate_adjacency_matrix(matrix,1)
    
    W_tmp = numpy.array(W)
    t1 = time.time()
    L_PYTHON = SpectralClusteringAlgorithm.calculate_laplacian(W_tmp, matrix, "PYTHON")[0]
    t2 = time.time()
    print 'Calculations took %0.3f s' % (t2-t1)
    
    W_tmp = numpy.array(W)
    t1 = time.time()
    L_NUMPY = SpectralClusteringAlgorithm.calculate_laplacian(W_tmp, matrix, "NUMPY")[0]
    t2 = time.time()
    print 'Calculations took %0.3f s' % (t2-t1)
    
    W_tmp = numpy.array(W)
    t1 = time.time()
    L_NUMPY_PURE = SpectralClusteringAlgorithm.calculate_laplacian(W_tmp, matrix, "NUMPY_PURE")[0]
    t2 = time.time()
    print 'Calculations took %0.3f s' % (t2-t1)
    
#     for i in range(len(L_PYTHON.flatten())):
#         print i
#         numpy.testing.assert_almost_equal(L_PYTHON.flatten()[i],L_NUMPY.flatten()[i],2)
    
    numpy.testing.assert_almost_equal(L_NUMPY, L_NUMPY_PURE, 10)
    numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY, 3)
    numpy.testing.assert_almost_equal(L_PYTHON, L_NUMPY_PURE, 3)
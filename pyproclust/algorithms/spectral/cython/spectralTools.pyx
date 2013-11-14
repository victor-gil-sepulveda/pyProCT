import pyproclust.algorithms.dbscan.cython.cythonDbscanTools as dbscanTools
import numpy
cimport numpy
import cython
import math

DOUBLE = numpy.double
ctypedef numpy.double_t DOUBLE_t


FLOAT = numpy.float
ctypedef numpy.float_t FLOAT_t


@cython.boundscheck(False)
def local_sigma_W_estimation(matrix):
    """
    Calculates local sigma estimation following Zelnik and Perona 
    
    @param matrix: The distance matrix for this dataset.
    
    @return: The adjacency matrix with the chosen sigma estimation.
    
    """
    cdef int N = matrix.row_length
    cdef int K = 7 # as stated in the paper
    cdef int i = 0
    cdef int j = 0
    cdef numpy.ndarray[FLOAT_t, ndim = 2] W_tmp = numpy.empty((N,N), dtype=FLOAT)
    cdef numpy.ndarray[FLOAT_t, ndim = 1] sigma = numpy.empty(N, dtype=FLOAT)
    cdef numpy.ndarray[FLOAT_t, ndim = 1] buffer = numpy.empty(N, dtype=FLOAT)
    
    # Finding local sigmas
    for i in range(N):
        sigma[i] = dbscanTools.kth_elements_distance(i, numpy.array([K]), buffer, matrix)[0]
    
    del buffer
    print "local sigmas found"
    # generating W
    for i in range(N):
        for j in range(i, N):
            W_tmp[i,j] = math.exp(-(matrix[i,j]*matrix[i,j])/(sigma[i]*sigma[j]))
            W_tmp[j,i] = W_tmp[i,j]
    #W_tmp = numpy.exp(-W_tmp)
    print "adjacency matrix created"
    return W_tmp, sigma
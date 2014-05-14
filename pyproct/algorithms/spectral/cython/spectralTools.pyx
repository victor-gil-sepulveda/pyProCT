import pyproct.algorithms.dbscan.cython.cythonDbscanTools as dbscanTools
import numpy
cimport numpy
import math
import scipy.linalg
import scipy.sparse.linalg
import cython
from pyRMSD.condensedMatrix import CondensedMatrix

DOUBLE = numpy.double
ctypedef numpy.double_t DOUBLE_t

FLOAT = numpy.float
ctypedef numpy.float_t FLOAT_t

MIN_SIGMA = 0.001

def order_by_eigenvalues(eigenvectors, eigenvalues):
    """
    Returns only a part.
    """
    # Order eigenvectors by eigenvalue (from lowest to biggest)
    idx = eigenvalues.real.argsort()
    return eigenvectors[:, idx]

@cython.boundscheck(False)
def calculate_degree_matrix(W):
    """
    Calculates... the degree matrix (as an array containing the diagonal).

    @param W: The adjacency matrix.

    @return: The degree matrix as an array.
    """
    # The degree matrix is the diagonal matrix with the degrees for each element.
    # The degree is the sum of the row of W (or column as it's symmetric) for each element.
    cdef int N = W.row_length
    cdef numpy.ndarray[FLOAT_t, ndim = 1] D = numpy.zeros(N, dtype=FLOAT)
    cdef double tmpdouble
    # Diagonal is 0 instead of 1, this means we can safely add it.
    for i in range(0,N):
        for j in range(0, N):
            tmpdouble = W[i,j]
            D[i] += tmpdouble
    return D

@cython.boundscheck(False)
def local_sigma_estimation(matrix):
    """
    Calculates local sigma estimation following Zelnik and Perona

    @param matrix: The distance matrix for this dataset.

    @return: The adjacency matrix with the squared mean of the sigma estimations.

    """
    cdef int N = matrix.row_length
    cdef int K = 7 # as stated in the paper
    cdef int i = 0
    cdef int j = 0

    cdef numpy.ndarray[FLOAT_t, ndim = 1] sigma = numpy.empty(N, dtype=FLOAT)
    cdef numpy.ndarray[FLOAT_t, ndim = 1] buffer = numpy.empty(N, dtype=FLOAT)

    # Finding local sigmas
    for i in range(N):
        sigma[i] = max(dbscanTools.kth_elements_distance(i, numpy.array([K]), buffer, matrix)[0],MIN_SIGMA)

    return sigma #numpy.mean(sigmas)**2

@cython.boundscheck(False)
def calculate_fully_connected_adjacency_matrix( matrix, sigma_sq):
    """
    Calculates the adjacency matrix (W) for the fully connected graph representation of this dataset.
    This matrix is modeled with a condensed matrix. Its diagonal should be 1 but due to the CondensedMatrix
    object specifications, it would be 0. It is of utmost importance to be aware of this when using it as part
    of other computations.

    @param matrix: The distance matrix for this dataset.

    @param sigma_sq: Is a cutoff parameter for edge creation.

    @return: The adjacency matrix.
    """
    W_data = numpy.exp(-((matrix.get_data()**2) / float(sigma_sq)))
    W = CondensedMatrix(W_data)
    return W

@cython.boundscheck(False)
def calculate_fully_connected_adjacency_matrix_with_sigma_estimation( matrix, sigmas):
    """

    @param matrix: The distance matrix for this dataset.

    @param sigmas:

    @return: The adjacency matrix.
    """
    cdef int N = matrix.row_length
    cdef numpy.ndarray[FLOAT_t, ndim = 1] W_data = numpy.empty((N*(N-1))/2, dtype=FLOAT)
    cdef int accum = 0

    for i in range(N-1):
        for j in range(i+1, N):
            W_data[accum] = math.exp(-(matrix[i,j]*matrix[i,j])/(sigmas[i]*sigmas[j]))
            accum = accum + 1
    W = CondensedMatrix(W_data)
    return W

def force_sparsity(W, sparsity = 0.5):
    """
    matrix is a CondensedMatrix
    """
    data = W.get_data()
    N = W.row_length
    # Find the weight that converts at least the 51% of the matrix in 0s
    # TODO: here we need an extra copy because we are directly working with the matrix memory, is this a bug
    # or a feature?
    threshold = sorted(data)[(((N*(N-1))/2) * sparsity)+1]

    for i in range(N-1):
        for j in range(i+1,N):
            if W[i,j] <= threshold:
                W[i,j] = 0.

@cython.boundscheck(False)
def calculateUnnormalizedLaplacian(W, D):
    """
    Implementation of the unnormalized laplacian calculation.
    Solves  L = D - W
    L is represented as a condensed matrix. As it needs to hold the diagonal, its number of elements
    must be N+1.

    @param W: Adjacency matrix. Its overwritten.
    @param D: Degree matrix.

    @return: Laplacian matrix (Symm. pos. definite).
    """

    cdef int N = W.row_length
    cdef numpy.ndarray[FLOAT_t, ndim = 2] L = numpy.zeros((N,N), dtype=FLOAT)

    for i in range(0,N):
        L[i][i] = D[i] - 1. # W[i][i] == 1

    for i in range(0,N-1):
        for j in range(i+1,N):
            tmpval = -W[i,j]
            L[i][j] = tmpval
            L[j][i] = tmpval
    return L

<<<<<<< HEAD
@cython.boundscheck(False)
=======
>>>>>>> 71b91d17cd806407651db13885ec3da4f14da1f5
def calculateUnnormalizedEigenvectors(L, max_number_of_clusters, is_sparse):
    if is_sparse:
        w, v  = scipy.sparse.linalg.eigsh(  A =L,
                                            k = max_number_of_clusters,
                                            which = 'SM',
                                            return_eigenvectors = True)
    else:
        w, v = scipy.linalg.eigh(a = L,
                                 type = 1,
                                 eigvals = (0, max_number_of_clusters-1),
                                 overwrite_a = True,
                                 check_finite = False)

    ordered_eigvectors = order_by_eigenvalues(v, w)
    return v


<<<<<<< HEAD
@cython.boundscheck(False)
=======

>>>>>>> 71b91d17cd806407651db13885ec3da4f14da1f5
def calculateNormalizedEigenvectors(L, D, max_number_of_clusters, is_sparse):
    if is_sparse:
        L_sparse = scipy.sparse.csc_matrix(L)
        w, v  = scipy.sparse.linalg.eigs(A = L,
                                         k = max_number_of_clusters,
                                         M = scipy.sparse.diags([D],[0],format = "csc"),
                                         which ='SM',
                                         return_eigenvectors = True)
    else:
        w, v = scipy.linalg.eigh(a=L,
                                 b= numpy.diag(D),
                                 eigvals = (0, max_number_of_clusters-1),
                                 overwrite_a = True,
                                 overwrite_b = True,
                                 check_finite = False)

    ordered_eigvectors = order_by_eigenvalues(v, w)
    return v

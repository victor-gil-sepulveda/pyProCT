from itertools import product
import pyproct.clustering.algorithms.spectral.cython.spectralTools as SpectralTools

def getClusterAndComplementary(i,all_clusters):
    """
    Returns the list representing the elements of cluster 'ith' and the
    list representing the sum of elements of all other clusters.
    """
    complementary_elements = []
    for j in range(len(all_clusters)):
        if j != i:
            complementary_elements.extend(all_clusters[j].all_elements)
    return(all_clusters[i].all_elements, complementary_elements)

cdef double W_partition (A1, A2, W):
    """
    W(A,B) = sum_{i pert A, j pert B} w_i_j
    """
    cdef double  w = 0.

    for i,j in product(A1,A2):
        w = w + W[i,j]
    return w

cdef double vol(A, D, condensed_matrix):
    """
    vol(A) = sum_{i pert A} d_i
    """
    cdef double vol_val = 0.
    cdef int i
    for i in  A:
        vol_val += D[i]
    return vol_val

cdef class CythonNCut(object):
    """
    Is a very simple separation measure.
    """
    def __init__(self):
        pass

    cpdef double evaluate(self, clustering, condensed_matrix) except *:
        """
        Evaluates:
        NCut = 1/2 sum_{i=1}^k W(A_i,A_i-complementary) / vol(A_i)
        """
        # Calculate similarity graph
        sigmas = SpectralTools.local_sigma_estimation(condensed_matrix)
        W  = SpectralTools.calculate_fully_connected_adjacency_matrix_with_sigma_estimation(condensed_matrix, sigmas)
        D = SpectralTools.calculate_degree_matrix(W)
        clusters = clustering.clusters

        cdef double ncut_val = 0
        cdef int i = 0
        for i in range(len(clusters)):
            A, Acomp = getClusterAndComplementary(i,clusters)
            ncut_val += W_partition(A, Acomp, condensed_matrix) / vol(A, D, condensed_matrix)
        return 0.5*ncut_val
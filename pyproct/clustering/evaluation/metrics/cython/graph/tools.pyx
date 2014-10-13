
from itertools import product
from pyproct.clustering.algorithms.spectral.cython import spectralTools
cimport cython

"""
Application of graph theory's metrics to clustering, where each connected component is
a cluster.
"""

def calc_adjacency_matrix(condensed_matrix):
    # Calculate adjacency matrix
    sigmas = spectralTools.local_sigma_estimation(condensed_matrix)
    W  = spectralTools.calculate_fully_connected_adjacency_matrix_with_sigma_estimation(condensed_matrix, sigmas)
    D = spectralTools.calculate_degree_matrix(W)
    return W, D

@cython.boundscheck(False)
def cut (A1, A2, adjacency_matrix):
    """
    W(A,B) = sum_{i pert A, j pert B} w_i_j in Luxburg's, cut in Shi & Malik.
    Sum of the weights of the edges one has to remove in order to separate two 
    components of the graph (it is 0 if components are not connected). 
    
    """
    cdef double  w = 0.

    for i,j in product(A1, A2):
        w = w + adjacency_matrix[i,j]
    return w

@cython.boundscheck(False)
def internal_vol(A, D, adjacency_matrix):
    """
    vol(A) = sum_{i pert A} d_i -> Luxburg's
    
    Clear graphical explanation in http://www.cis.upenn.edu/~jshi/GraphTutorial/Tutorial-ImageSegmentationGraph-cut1-Shi.pdf
    
    The volume is the sum of all the weights of the edges of a partition (including those that connect with other partitions), but 
    each edge must be counted only once.
    
    Internal edges are counted twice, so its contribution must be multiplied by 0.5.
    
    This function calculates only the internal volume.  
    """
    
    cdef double internal_vol = 0.
    cdef int i
    
    for i in A:
        internal_vol += D[i]
    
    return 0.5 * internal_vol 

@cython.boundscheck(False)
def get_cluster_and_complementary (i, all_clusters):
    """
    Returns a list with the elements of 'ith' cluster and anoher
    list with the elements of all other clusters.
    :param i: The index of the cluster we want to use in 'all_clusters'.
    :param all_clusters: A list with all the clusters.
    
    :return: A list with ith cluster elements and its complementary.
    """
    complementary_elements = []
    
    for j in range(len(all_clusters)):
        if j != i:
            complementary_elements.extend(all_clusters[j].all_elements)
            
    return all_clusters[i].all_elements, complementary_elements
            
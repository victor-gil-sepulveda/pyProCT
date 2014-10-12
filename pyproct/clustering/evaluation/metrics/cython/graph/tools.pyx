
from itertools import product
from pyproct.clustering.algorithms.spectral.cython import spectralTools

"""
Application of graph theory's metrics to clustering, where each connected component is
a cluster.
"""

cdef double cut (A1, A2, adjacency_matrix):
    """
    W(A,B) = sum_{i pert A, j pert B} w_i_j in Luxburg's, cut in Shi & Malik.
    Sum of the weights of the edges one has to remove in order to separate two 
    components of the graph (it is 0 if components are not connected). 
    
    """
    cdef double  w = 0.

    for i,j in product(A1, A2):
        w = w + adjacency_matrix[i,j]
    return w

cdef double internal_vol(A, D, adjacency_matrix):
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
            
            
            
            

def d(i, adjacency_condensed_matrix):
    """
    Degree of a vertex:
    d_i = sum_{j=1}^n w_ij
    """
    d_val = 0
    n = condensed_matrix.row_length
    for j in range(n):
        d_val += adjacency_condensed_matrix[i,j]
    return d_val


def vol(A, adjacency_condensed_matrix):
    """
    vol(A) = sum_{i pert A} d_i
      where the degree matrix is the degree matrix of the unconnected graph.
    """
    vol_val = 0
    for i in  A:
        vol_val += d(i, adjacency_condensed_matrix)   
    return vol_val

def W_partition (A1, A2, W):
    """
    W(A,B) = sum_{i pert A, j pert B} w_i_j
    """
    w = 0.

    for i,j in product(A1,A2):
        w = w + W[i,j]
    return w

def W(A1, A2, condensed_matrix):
    """
    W(A,B) = sum_{i pert A, j pert B} w_i_j
    """
    w = 0;
    
    for indices in product(A1,A2):
        w = w + condensed_matrix[indices]
    return w

def vol2(A, D, condensed_matrix):
    """
    vol2(A) = sum_{i pert A} d_i
    """
    vol_val = 0
    for i in  A:
        vol_val += D[i]
    return vol_val

def all_cut(clustering,condensed_matrix):
    """
    Cut measure of the connected component A vs the complementary set.
    """
    clusters = clustering.clusters
    cut_val = 0
    for i in range(len(clusters)):
        A,Acomp = getClusterAndComplementary(i,clusters)
        cut_val += W(A,Acomp,condensed_matrix)
    return 0.5*cut_val

def single_cut(A,Acomp,condensed_matrix):
    return 0.5 * W(A,Acomp,condensed_matrix)
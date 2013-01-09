from itertools import product

def getClusterAndComplementary(i,all_clusters):
    """
    Returns the list representing the elements of cluster 'ith' and the 
    list representing the sum of elements of all other clusters. 
    """
    complementary = []
    for j in range(len(all_clusters)):
        if j != i:
            complementary.extend(all_clusters[j].all_elements)
    return(all_clusters[i].all_elements,complementary)

cdef double W(A1,A2,condensed_matrix):
    """
    W(A,B) = sum_{i pert A, j pert B} w_i_j
    """
    cdef double  w = 0.
    
    for indices in product(A1,A2):
        w = w + condensed_matrix[indices]
    return w

cdef double d( int i, condensed_matrix):
    """
    Degree of a vertex:
    d_i = sum_{j=1}^n w_ij
    """
    cdef double d_val = 0.
    cdef int n = condensed_matrix.row_length
    cdef int j
    for j in range(n):
        d_val += condensed_matrix[i,j]
    return d_val

cdef double vol(A,condensed_matrix):
    """
    vol(A) = sum_{i pert A} d_i
    """
    cdef double vol_val = 0.
    cdef int i
    for i in  A:
        vol_val += d(i,condensed_matrix)   
    return vol_val


cdef class CythonNCut(object):
    """
    Is a very simple separation measure.
    """
    def __init__(self):
        pass
    
    cpdef double evaluate(self,clustering,condensed_matrix) except *:
        """
        Evaluates:
        NCut = 1/2 sum_{i=1}^k W(A_i,A_i-complementary) / vol(A_i)
        """
        clusters = clustering.clusters
        cdef double ncut_val = 0
        cdef int i = 0
        for i in range(len(clusters)):
            A, Acomp = getClusterAndComplementary(i,clusters)
            ncut_val += W(A, Acomp, condensed_matrix) / vol(A, condensed_matrix)
        return 0.5*ncut_val
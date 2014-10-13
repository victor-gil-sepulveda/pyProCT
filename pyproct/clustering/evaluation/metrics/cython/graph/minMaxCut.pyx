from pyproct.clustering.evaluation.metrics.cython.graph.tools import get_cluster_and_complementary, cut,\
                                            calc_adjacency_matrix, internal_vol

cdef class MinMaxCut(object):
    """
    Within cluster similarity.
    """
    def __init__(self):
        pass
    
    cpdef double evaluate(self, clustering, condensed_matrix) except *:
        """
        Evaluates:
        MinMaxCut = sum_{i=1}^k cut(A_i,A_i-complementary) / W(A_i,A_i)
        """
        clusters = clustering.clusters
        W, D = calc_adjacency_matrix(condensed_matrix)
        
        cdef double mmcut_val = 0
        
        for i in range(len(clusters)):
            A, Acomp = get_cluster_and_complementary(i, clusters)
            mmcut_val += cut(A, Acomp, W) / internal_vol(A, D, condensed_matrix)
        
        return 0.5*mmcut_val
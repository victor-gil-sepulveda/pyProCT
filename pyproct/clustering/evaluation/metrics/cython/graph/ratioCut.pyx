from pyproct.clustering.evaluation.metrics.cython.graph.tools import get_cluster_and_complementary, cut,\
                                            calc_adjacency_matrix

cdef class RatioCut(object):
    """
    Is another simple separation measure.
    """
    def __init__(self):
        pass
    
    cpdef double evaluate(self, clustering, condensed_matrix) except *:
        """
        Evaluates:
        RatioCut = sum_{i=1}^k cut(A_i, A_i_complementary) / card(A_i)
        """
        clusters = clustering.clusters
        W,D = calc_adjacency_matrix(condensed_matrix)
        
        cdef double ratiocut_val = 0
        
        for i in range(len(clusters)):
            A, Acomp = get_cluster_and_complementary(i, clusters)
            ratiocut_val += cut(A, Acomp, W) / len(A)
        return ratiocut_val
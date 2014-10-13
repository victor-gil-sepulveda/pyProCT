from itertools import product
from pyproct.clustering.evaluation.metrics.cython.graph.tools import calc_adjacency_matrix, cut,\
                                                internal_vol , get_cluster_and_complementary


cdef class NCut(object):
    """
    Separation measure. Calculates the average percent of edges (weights) that separate the partitions with
    the weights of the  partitions themselves.
    
    Formal descriptions in:
        Shi and Malik 2000
        Luxburg 2007
    """
    def __init__(self):
        pass
    
    cpdef double evaluate(self, clustering, condensed_matrix) except *:
        """
        Evaluates:
        NCut = 1/2 sum_{i=1}^k W(A_i,A_i_complementary) / vol(A_i) -> Luxburg's
                        
        :param clustering: The clustering to compute the measure.
        :param adjacency_condensed_matrix: The adjacency matrix of the graph representation of the 
        clustering (elements act as nodes and distances as weights, a kernel is used to express
        connectivity).
        
        :return: The value of the measure (double).
        """
        # Calculate adjacency matrix
        clusters = clustering.clusters
        W,D = calc_adjacency_matrix(condensed_matrix)

        cdef double ncut_val = 0
        cdef int i = 0
        cdef double cut_value, internal_volume
        
        for i in range(len(clusters)):
            A, Acomp = get_cluster_and_complementary(i, clusters)
            cut_value =  cut(A, Acomp, W)
            internal_volume = internal_vol(A, D, W)
            volume = internal_volume + cut_value
            ncut_val += cut_value / volume

        # Divided by K in http://www.cis.upenn.edu/~jshi/GraphTutorial/Tutorial-ImageSegmentationGraph-cut1-Shi.pdf
        return ncut_val / len(clusters)    
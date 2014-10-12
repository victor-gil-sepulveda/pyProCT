from pyproct.clustering.evaluation.metrics.cython.graph.nCut import get_cluster_and_complementary



class RatioCut(object):
    """
    Is another simple separation measure.
    """
    def __init__(self):
        pass
    
    def evaluate(self, clustering, condensed_matrix):
        """
        Evaluates:
        RatioCut = sum_{i=1}^k cut(A_i,A_i-complementary) / card(A_i)
        """
        clusters = clustering.clusters
        ratiocut_val = 0
        for i in range(len(clusters)):
            A, Acomp = get_cluster_and_complementary(i,clusters)
            
            ratiocut_val += W(A,Acomp,condensed_matrix) / len(A)
        return 0.5*ratiocut_val
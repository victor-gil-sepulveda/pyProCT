class MinMaxCut(object):
    """
    Within cluster similarity.
    """
    def __init__(self):
        pass
    
    def evaluate(self, clustering, condensed_matrix):
        """
        Evaluates:
        MinMaxCut = sum_{i=1}^k cut(A_i,A_i-complementary) / W(A_i,A_i)
        """
        clusters = clustering.clusters
        mmcut_val = 0
        for i in range(len(clusters)):
            A, Acomp = getClusterAndComplementary(i, clusters)
            mmcut_val += single_cut(A,Acomp,condensed_matrix) / W(A,A,condensed_matrix)
        return 0.5*mmcut_val
"""
Created on 13/08/2012

@author: victor
"""

from itertools import product

"""
Application of graph theory's metrics to clustering, where each connected component is
a cluster.
"""
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

def W(A1,A2,condensed_matrix):
    """
    W(A,B) = sum_{i pert A, j pert B} w_i_j
    """
    w = 0;
    
    for indices in product(A1,A2):
        w = w + condensed_matrix[indices]
    return w

def d(i,condensed_matrix):
    """
    Degree of a vertex:
    d_i = sum_{j=1}^n w_ij
    """
    d_val = 0
    n = condensed_matrix.row_length
    for j in range(n):
        d_val += condensed_matrix[i,j]
    return d_val

def vol(A,condensed_matrix):
    """
    vol(A) = sum_{i pert A} d_i
    """
    vol_val = 0
    for i in  A:
        vol_val += d(i,condensed_matrix)   
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
    
class MinMaxCut(object):
    """
    Within cluster similarity.
    """
    def __init__(self):
        pass
    
    def evaluate(self,clustering,condensed_matrix):
        """
        Evaluates:
        MinMaxCut = sum_{i=1}^k cut(A_i,A_i-complementary) / W(A_i,A_i)
        """
        clusters = clustering.clusters
        mmcut_val = 0
        for i in range(len(clusters)):
            A,Acomp = getClusterAndComplementary(i,clusters)
            mmcut_val += single_cut(A,Acomp,condensed_matrix) / W(A,A,condensed_matrix)
        return 0.5*mmcut_val

class NCut(object):
    """
    Is a very simple separation measure.
    """
    def __init__(self):
        pass
    
    def evaluate(self,clustering,condensed_matrix):
        """
        Evaluates:
        NCut = 1/2 sum_{i=1}^k W(A_i,A_i-complementary) / vol(A_i)
        """
        clusters = clustering.clusters
        ncut_val = 0
        for i in range(len(clusters)):
            A,Acomp = getClusterAndComplementary(i,clusters)
            ncut_val += W(A,Acomp,condensed_matrix) / vol(A,condensed_matrix)
        return 0.5*ncut_val

class RatioCut(object):
    """
    Is another simple separation measure.
    """
    def __init__(self):
        pass
    
    def evaluate(self,clustering,condensed_matrix):
        """
        Evaluates:
        RatioCut = sum_{i=1}^k cut(A_i,A_i-complementary) / card(A_i)
        """
        clusters = clustering.clusters
        ratiocut_val = 0
        for i in range(len(clusters)):
            A,Acomp = getClusterAndComplementary(i,clusters)
            ratiocut_val += W(A,Acomp,condensed_matrix) / len(A)
        return 0.5*ratiocut_val
    

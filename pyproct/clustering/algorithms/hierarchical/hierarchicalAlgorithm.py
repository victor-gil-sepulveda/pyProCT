"""
Created on 19/03/2012

@author: victor
"""
import scipy.cluster.hierarchy as hcluster
import  fastcluster as hcluster_fast
from pyproct.clustering.cluster import gen_clusters_from_class_list
from pyproct.clustering.clustering import Clustering

class HierarchicalClusteringAlgorithm(object):
    """
    TODOC
    """
    def __init__(self,condensed_matrix, **kwargs):
        """
        Constructor.
        """
        self.condensed_matrix = condensed_matrix
        self.hie_mat = None

    def perform_clustering(self, kwargs):
        """
        Performs the hierarchical clustering step and the clustering step. If the hierarchical
        matrix is given, then it just calculates the clusters for a given cutoff. If we call the algorithm
        a second time it will use the last matrix.
        """
        """
        Gets a condensed matrix and calculates the clustering. One can use
        diverse methodologies to do this clustering...
        With preserve_input=False the matrix is destroyed while clustering, ut it saves
        memory.
        The metric is not needed in this case,as we are giving the function the calculated
        matrix.
        The method is the method used to determine distances when fusing clusters. methods are described in:
        http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
        """
        try:
            cutoff = kwargs["cutoff"]
        except KeyError:
            cutoff = None

        try:
            hie_mat = kwargs["hie_mat"]
        except KeyError:
            hie_mat = None

        try:
            method = kwargs["method"]
        except KeyError:
            method = 'complete'

        if hie_mat != None:
            self.hie_mat = hie_mat
#            print "[HIERARCHICAL] Matrix provided."
        else:
            if self.hie_mat == None:
                #self.hie_mat = fast_hcluster.linkage(condensed_matrix, method='centroid', metric='euclidean', preserve_input=False)
#                print "[HIERARCHICAL] Calculating Matrix"
                #self.hie_mat = fastclust.linkage(self.condensed_matrix.get_data(), method = method)
                self.hie_mat = hcluster_fast.linkage(self.condensed_matrix.get_data(), method = method)
#            else:
#                print "[HIERARCHICAL] Matrix was already stored"

        algorithm_details = "Hierarchical with "+method+" method (cutoff = " +str(cutoff)+")"

        if cutoff != None:
            # Then apply the cutoff, this doesn't work much as expected
#            print "[HIERARCHICAL] getting clustering."+algorithm_details
            group_list = hcluster.fcluster(self.hie_mat,cutoff)
#            print "[HIERARCHICAL] Clustering done."+algorithm_details
            # Then let's generate the clusters
            clusters = gen_clusters_from_class_list(group_list)
            return Clustering(clusters,details = algorithm_details)
        else:
            return None

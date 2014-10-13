"""
Created on 06/06/2013

@author: victor
"""
import numpy
from pyproct.clustering.evaluation.metrics.common import get_inter_cluster_distances,\
    get_intra_cluster_distances

class DunnCalculator(object):
    """
    Calculates Dunn index as presented in Krzystof Kryszczuk and Paul Hurley's paper . 
    """
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        dmin = self.min_intracluster_distances(clustering, matrix)
        dmax = self.max_intercluster_distance(clustering, matrix)
        dunn_index = dmin / dmax
        return dunn_index
    
    @classmethod
    def min_intracluster_distances(cls, clustering, matrix):
        """
        Calculates d_min, the minimum internal distance.
        @param clustering: The clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: d_min's value
        """
        return numpy.min([numpy.min(get_intra_cluster_distances(c, matrix)) for c in clustering.clusters])
    
    @classmethod
    def max_intercluster_distance(cls, clustering, matrix):
        """
        Calculates d_max, the maximum inter clusters.
        @param clustering: The clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: d_max' value
        """
        distances = []
        for i in range(len(clustering.clusters)-1):
            for j in range(i+1,len(clustering.clusters)):
                distances.extend(get_inter_cluster_distances( i, j, clustering.clusters, matrix))
        return numpy.max(distances)
    
    
    
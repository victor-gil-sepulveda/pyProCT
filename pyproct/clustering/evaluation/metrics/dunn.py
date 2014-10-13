"""
Created on 06/06/2013

@author: victor
"""
import numpy
from pyproct.clustering.evaluation.metrics.common import get_inter_cluster_distances,\
    get_intra_cluster_distances
from pyproct.tools.exceptions import SingularClusterException

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
        distances = []
        for c in clustering.clusters:
            try:
                distances.append(numpy.min(get_intra_cluster_distances(c, matrix)))
            except SingularClusterException:
                # If we work with a singular cluster, we add 0s so that no min function
                # fails. The convention for the distance of a cluster with only one element
                # will be 0 in this case.
                distances.append(0)
        return numpy.min(distances)
    
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
    
    
    
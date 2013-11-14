'''
Created on 06/06/2013

@author: victor
'''
import numpy
from pyproclust.clustering.metrics.common import get_distances_of_elements_to,\
    update_medoids
from pyproclust.clustering.cluster import Cluster

class CompactnessCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        """
        mean is approximated to medoid 
        """
        update_medoids(clustering, matrix)
        
        # Create a cluster that contains all elements to calculate variances, we'll check even 
        # non clustered elements for the sake of simplicity
        global_cluster = Cluster(None, clustering.get_all_clustered_elements())
        global_cluster.prototype = global_cluster.calculate_medoid(matrix)
        global_variance = numpy.var(get_distances_of_elements_to(global_cluster.prototype, 
                                                                 global_cluster.all_elements, 
                                                                 matrix))
        variances = [self.cluster_variance(cluster,matrix) for cluster in clustering.clusters]
        
        sum_ci = numpy.sum(variances)
        
        Cmp = sum_ci / (len(clustering.clusters)*global_variance)
        
        return Cmp
    
    @classmethod
    def cluster_variance(cls, cluster, matrix):
        """
        precondition, cluster medoid (prototype) it's alread
        """
        return numpy.var(get_distances_of_elements_to(cluster.prototype, cluster.all_elements, matrix))

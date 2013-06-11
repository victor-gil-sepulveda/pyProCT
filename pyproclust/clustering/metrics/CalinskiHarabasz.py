'''
Created on 11/06/2013

@author: victor
'''
import numpy
from pyproclust.clustering.metrics.common import update_medoids,  get_intra_cluster_distances

class CalinskiHarabaszCalculator(object):
    
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        update_medoids(clustering, matrix)
        # We'll follow the paper expanded formula
        
        # We need the general mean
        global_mean_distance = numpy.mean(matrix.get_data())
        print 
        print "d", global_mean_distance
        # A_k calculation
        k = len(clustering.clusters)
        n = matrix.row_length
        print [self.ch_cluster_term(cluster, global_mean_distance, matrix) for cluster in clustering.clusters]
        A_k = (1./( n - k))*numpy.array([self.ch_cluster_term(cluster, global_mean_distance, matrix) for cluster in clustering.clusters]).sum();
        print A_k
        VRC = (global_mean_distance + ((n -k) / (k-1))*A_k) / float((global_mean_distance - A_k))
        
        return VRC
    @classmethod
    def ch_cluster_term(cls, cluster, global_mean_distance, matrix):
        # Calculate cluster mean distance
        n = len(cluster.all_elements)
        cluster_mean_distance = numpy.mean(numpy.array(get_intra_cluster_distances( cluster, matrix)))
        print
        print "cluster_mean_dist", cluster_mean_distance, get_intra_cluster_distances(cluster, matrix)
        return (n-1) * (global_mean_distance - cluster_mean_distance)
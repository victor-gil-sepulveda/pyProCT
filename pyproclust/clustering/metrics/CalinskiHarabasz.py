'''
Created on 11/06/2013

@author: victor
'''
import numpy
from pyproclust.clustering.metrics.common import update_medoids,  get_intra_cluster_distances

def mean(array):
    if len(array)==0:
        return 0.
    else:
        return numpy.mean(array)

class CalinskiHarabaszCalculator(object):
    
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        # Cluster prototypes (medoids here) must be updated
        update_medoids(clustering, matrix)
        
        # We'll follow the paper expanded formula
        # First we need the general mean of the squared distances
        D = mean(matrix.get_data()**2)
        # A_k calculation
        k = len(clustering.clusters)
        n = matrix.row_length
        print "WGSS", self.WGSS(clustering.clusters, matrix)
        print "BGSS", self.BGSS(clustering, D,  matrix)
        print (self.BGSS(clustering, D,  matrix)/self.WGSS(clustering.clusters, matrix))*(float(n-k)/(k-1))
        A_k = CalinskiHarabaszCalculator.A_k(clustering, D, matrix)
        VRC = (D + (float(n -k) / (k-1))*A_k) / float((D - A_k))
        return VRC
    
    
    @classmethod
    def WGSS(cls, clusters, matrix):
        wgss = 0
        for c in clusters:
            n = len(c.all_elements)
            d = mean(numpy.array(get_intra_cluster_distances( c, matrix))**2) 
            wgss += (n-1)*d
        return wgss*0.5
    
    @classmethod
    def BGSS(cls, clustering, D, matrix):
        n = clustering.total_number_of_elements
        k = len(clustering.clusters)
        bgss = (k-1)*D + (n-k)*cls.A_k(clustering, D, matrix)
        return bgss*0.5
    
    @classmethod
    def A_k(cls, clustering, D, matrix):
        n = clustering.total_number_of_elements
        k = len(clustering.clusters)
        return (1./( n - k))*numpy.array([cls.ch_cluster_term(cluster, D, matrix) for cluster in clustering.clusters]).sum();
    
    @classmethod
    def ch_cluster_term(cls, cluster, global_mean_distance, matrix):
        """
        Calculates one of the formula terms (ng-1)(D-d_g)
        @param cluster: The cluster to use in calculation.
        @param global_mean_distance: 'D'. Is the mean of the n*(n-1)/2 distances of all the elements.
        @
        """
        # Calculate cluster mean distance
        n = len(cluster.all_elements)
        cluster_mean_distance = mean(numpy.array(get_intra_cluster_distances( cluster, matrix))**2)
        return (n-1) * (global_mean_distance - cluster_mean_distance)
"""
Created on 11/06/2013

@author: victor
"""
import numpy
from pyproct.clustering.evaluation.metrics.common import update_medoids,  get_intra_cluster_distances

def mean(array):
    if len(array)==0:
        return 0.
    else:
        return numpy.mean(array)

class CalinskiHarabaszCalculator(object):
    """
    Implementation of the Calinski-Harabasz (VRC) index, as described in T. Calinski and J.Harabasz paper in
    'Communications in Statistics 3(1), 1-27 (1974)'.
    Contains the necessary functions to calculate the index in the two forms presented.
    """
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        """
        Calculates the index value for a clustering.
        @param clustering: The clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: The calculated Calinski-Harabasz (VRC) index.
        """
        # Cluster prototypes (medoids here) must be updated
        update_medoids(clustering, matrix)
        
        # We'll follow the paper expanded formula
        # First we need the global mean of the squared distances
        D = mean(matrix.get_data()**2)
        # A_k calculation
        k = len(clustering.clusters)
        n = matrix.row_length
#         print "WGSS", self.WGSS(clustering.clusters, matrix)
#         print "BGSS", self.BGSS(clustering, D,  matrix)
#         print (self.BGSS(clustering, D,  matrix)/self.WGSS(clustering.clusters, matrix))*(float(n-k)/(k-1))
        A_k = CalinskiHarabaszCalculator.A_k(clustering, D, matrix)
        VRC = (D + (float(n -k) / (k-1))*A_k) / float((D - A_k))
        return VRC
    
    
    @classmethod
    def WGSS(cls, clusters, matrix):
        """
        C-H description of the "Within group sum of squares".
        @param clusters: An array with all clusters description (usually Clustering.clusters)
        @param matrix: The condensed matrix containing all distances.
        @return: The value of WGSS.
        """
        wgss = 0
        for c in clusters:
            n = len(c.all_elements)
            d = mean(numpy.array(get_intra_cluster_distances( c, matrix))**2) 
            wgss += (n-1)*d
        return wgss*0.5
    
    @classmethod
    def BGSS(cls, clustering, D, matrix):
        """
        C-H description of the "Between group sum of squares".
        @param clustering: The clustering being checked.
        @param D: Mean distance of the sum of all squared distances present in the matrix.
        @param matrix: The condensed matrix containing all distances.
        @return: The value of BGSS.
        """
        n = clustering.total_number_of_elements
        k = len(clustering.clusters)
        bgss = (k-1)*D + (n-k)*cls.A_k(clustering, D, matrix)
        return bgss*0.5
    
    @classmethod
    def A_k(cls, clustering, D, matrix):
        """
        Calculates the A_k term.
        @param clustering: The clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: The A_k term value.
        """
        n = clustering.total_number_of_elements
        k = len(clustering.clusters)
        return (1./( n - k))*numpy.array([cls.ch_cluster_term(cluster, D, matrix) for cluster in clustering.clusters]).sum();
    
    @classmethod
    def ch_cluster_term(cls, cluster, global_mean_distance, matrix):
        """
        Calculates one of the formula terms (ng-1)(D-d_g)
        @param cluster: The cluster to use in calculation.
        @param global_mean_distance: 'D'. Is the mean of the n*(n-1)/2 distances of all the elements.
        @param matrix: The condensed matrix containing all distances.
        @return: Calculated term.
        """
        # Calculate cluster mean distance
        n = len(cluster.all_elements)
        cluster_mean_distance = mean(numpy.array(get_intra_cluster_distances( cluster, matrix))**2)
        return (n-1) * (global_mean_distance - cluster_mean_distance)
    
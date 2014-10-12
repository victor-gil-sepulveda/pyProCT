"""
Created on May 30, 2013

@author: victor
"""
import numpy
from pyproct.clustering.evaluation.metrics.common import update_medoids,\
    get_distances_of_elements_to

class DaviesBouldinCalculator(object):
    """
    Is an estimator of cluster validity. The lower the better.
    """
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        """
        Calculates the index value for a clustering.
        @param clustering: The clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: The calculated Davies-Bouldin estimator value.
        """
        # Medoids will be the new prototypes
        update_medoids(clustering, matrix)
        
        # We calculate the average distances to the medoid for each of the clusters
        average_distances = self.calc_average_distances(clustering, matrix)
        
        # Then we can do the actual calculation
        db_index = self.db_index_calculation( average_distances, clustering.clusters, matrix)
        
        return db_index
    
    @classmethod
    def calc_average_distances(cls, clustering, matrix):
        return [cls.calculate_average_distance_from_prototype(cluster, matrix) for cluster in clustering.clusters]
    
    @classmethod
    def db_index_calculation(cls, d, clusters, matrix):
        """
        Calculates the summation of the terms for all clusters vs all the other clusters of the clustering.
        @param d: Array containing the mean distances to the centroid (medoid) for all clusters (with index
        correspondence).
        @param clusters: An array containing all the clusters of the clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: The index value.
        """
        db_index = 0
        k = len(clusters)
        for i in range(k-1):
            db_index += numpy.max(cls.db_terms_for_cluster(i, d, clusters, matrix))
        return db_index/k
        
    @classmethod
    def db_terms_for_cluster(cls, i, d, clusters, matrix):
        """
        Calculates the summation of the terms for a cluster vs all the other clusters of the clustering.
        @param i: Is the index of the first cluster (in 'clusters') to calculate the term.
        @param d: Array containing the mean distances to the centroid (medoid) for all clusters (with index
        correspondence).
        @param clusters: An array containing all the clusters of the clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: The value for each term (i vs j[0..n] j!=i).
        """
        terms = []
        for j in range(i+1, len(clusters)):
            if i != j:
                terms.append(cls.d_b_term_calculation(i,j,d,clusters, matrix))
        return terms
    
    @classmethod
    def d_b_term_calculation(cls, i, j, d, clusters, matrix):
        """
        Calculates one of the terms as defined in the paper (di+dj) / d(ci,cj)
        @param i: Is the index of the first cluster (in 'clusters') to calculate the term.
        @param j: The index of the second cluster.
        @param clusters: An array containing all the clusters of the clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: The value of the term.
        """
        c_i, c_j = clusters[i].prototype, clusters[j].prototype
        return (d[i]+d[j]) / matrix[c_i, c_j]
    
    @classmethod    
    def calculate_average_distance_from_prototype(cls, cluster, matrix):
        """
        Returns the average distance of the elements of a cluster with its medoid.
        @param cluster: The cluster from which we want to calculate this distance.
        @param matrix: The condensed matrix containing all distances.
        @return: The calculated value.
        """
        proto = cluster.prototype
        elements_copy = list(cluster.all_elements)
        elements_copy.remove(proto)
        distances = get_distances_of_elements_to(proto, elements_copy, matrix)
        if distances == []:
            return 0.
        else:
            return numpy.mean(distances)   
    
        
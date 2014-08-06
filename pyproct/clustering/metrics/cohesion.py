"""
Created on 09/01/2013

@author: victor
"""

class CohesionCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        total_cohesion = 0
        for cluster in clustering.clusters:
            total_cohesion = total_cohesion + self.evaluate_cluster(cluster, matrix)
        return total_cohesion
    
    def evaluate_cluster(self, cluster, condensed_distance_matrix):
        """
        Returns the cohesion value of a cluster. The condensed matrix given as 
        parameter stores the distances of the elements of the dataset used to extract
        the cluster.
        
        The definition of cohesion would be weight*2*cohesion in the case we follow
        the exact formula in the book []. As we are going to do comparisons, the x2 global
        multiplication doesn't affect.
        
        Cohesion of a cluster of 1 element should be infinite instead of 0...
        """
        size = cluster.get_size()
        weight = 1. / size
        cohesion = 0.
        for i in range(size-1):
            for j in range(i+1,size):
                cohesion = cohesion + condensed_distance_matrix[cluster[i],cluster[j]]
        return weight*cohesion
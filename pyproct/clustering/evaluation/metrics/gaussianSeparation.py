"""
Created on 06/06/2013

@author: victor
"""
from pyproct.clustering.evaluation.metrics.common import update_medoids,\
    get_inter_cluster_prototype_distances
import numpy

class GaussianSeparationCalculator(object):
    """
    Cluster separation index calculation, as defined in Wu and H. Xiong 2002.
    """
    TWO_SIGMA_SQUARED = 0.25
    
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        """
        Calculates the index.
        @param clustering: The clustering being checked.
        @param matrix: The condensed matrix containing all distances.
        @return: The calculated value for the index.
        """
        update_medoids(clustering, matrix)
        
        C = len(clustering.clusters)
        constant = 2. / (C * (C-1)) # x2 as we use only half of the distances  
        
        Sep = constant * self.exponential_list_generation(clustering, matrix).sum()
        
        return Sep
    
    @classmethod
    def exponential_list_generation(cls, clustering, matrix):
        # May not work if prototypes are not defined :S
        proto_distances = numpy.array(get_inter_cluster_prototype_distances(clustering.clusters, matrix))
        return numpy.exp(((proto_distances**2)/ -cls.TWO_SIGMA_SQUARED))
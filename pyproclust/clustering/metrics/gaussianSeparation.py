'''
Created on 06/06/2013

@author: victor
'''
from pyproclust.clustering.metrics.common import update_medoids
import numpy

class GaussianSeparationCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        """
        """
        TWO_SIGMA_SQUARED = 0.25
        
        constant = 1. / (len(clustering.clusters)* (len(clustering.clusters)-1))
        
        update_medoids(clustering, matrix)
        
        Sep = constant * numpy.exp(((self.inter_cluster_prototype_distances()**2)/-TWO_SIGMA_SQUARED)).sum()
        
        return Sep
'''
Created on 06/06/2013

@author: victor
'''
class GaussianSeparationCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        """
        """
        TWO_SIGMA_SQUARED = 0.25
        
        constant = 1. / (len(clustering.clusters)* (len(clustering.clusters)-1))
        
        Sep = constant * self.inter_cluster_gaussian_separation()
        
    @classmethod
    def inter_cluster_gaussian_separation(cls):
        pass
    
    @classmethod
    def inter_cluster_pairwise_gaussian_separation(cls, i, j, clusters, matrix):
        pass
        
    
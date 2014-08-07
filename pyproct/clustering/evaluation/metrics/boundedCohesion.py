"""
Created on 09/01/2013

@author: victor
"""
import numpy
from pyproct.tools.matrixTools import get_submatrix

class MirrorCohesionCalculator(object):
    """
    Bounded in [0,1] where 0 means the perfect cohesion and 1 the worst. Doesn't use the prototype.
    """
    def __init__(self):
        pass
    
    def evaluate(self, clustering, condensed_distance_matrix):
        """
        Returns the cohesion value of a cluster. The weight will be the number of elements
        of each cluster. 
        
        The maximum value will be (in the case we have only one cluster) the sum of
        all the pair distances. We can avoid the 2x factor, increasing the performance.
        """
        clustered_elements = sorted(clustering.get_all_clustered_elements())
        number_of_clustered_elements = len(clustered_elements)
        if number_of_clustered_elements > 0:
            distances = get_submatrix(condensed_distance_matrix, clustered_elements)
            max_cohesion = numpy.sum(distances.get_data())/number_of_clustered_elements
            
            total_cohesion = 0
            for c in clustering.clusters:
                size = c.get_size()
                weight = 1. / size
                cohesion = 0.
                for i in range(size-1):
                    for j in range(i+1, size):
                        cohesion = cohesion + condensed_distance_matrix[c[i],c[j]]
                total_cohesion +=  weight*cohesion
            return 1 - total_cohesion / max_cohesion
        else:
            return 0.

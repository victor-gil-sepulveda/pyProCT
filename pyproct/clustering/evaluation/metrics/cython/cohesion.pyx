# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
"""
Created on 12/07/2012

@author: victor
"""

import numpy
cimport numpy
cimport cython
numpy.import_array()
from pyproct.tools.matrixTools import get_submatrix

cdef class CohesionCalculator(object):
    def __init__(self):
        pass
    
    @cython.boundscheck(False)
    cpdef double evaluate(self, clustering, condensed_distance_matrix)  except *:
        cdef double weight,cohesion,max_cohesion
        cdef int size, i, j, total_number_of_elements
        
        clustered_elements = sorted(clustering.get_all_clustered_elements())
        number_of_clustered_elements = len(clustered_elements)
        if number_of_clustered_elements > 0:
            distances = get_submatrix(condensed_distance_matrix, clustered_elements)
            max_cohesion = numpy.sum(distances.get_data())/number_of_clustered_elements
            del distances
            
            total_cohesion = 0
            for c in clustering.clusters:
                size = c.get_size()
                weight = 1. / size
                cohesion = 0.
                for i in range(size-1):
                    for j in range(i+1,size):
                        cohesion = cohesion + condensed_distance_matrix[c[i],c[j]]
                total_cohesion +=  weight*cohesion
            return 1 - (total_cohesion / max_cohesion)
        else:
            return 0.
    
    @cython.boundscheck(False)
    cpdef double  evaluate_cluster(self, cluster, condensed_distance_matrix):
        """
        Returns the cohesion value of a cluster. The condensed matrix given as 
        parameter stores the distances of the elements of the dataset used to extract
        the cluster.
        
        The definition of cohesion would be weight*2*cohesion in the case we follow
        the exact formula in the book []. As we are going to do comparisons, the x2 global
        multiplication doesn't affect.
        
        Cohesion of a cluster of 1 element should be infinite instead of 0...
        """
        cdef int size = cluster.get_size()
        cdef double weight = 1. / size
        cdef int i = 0
        cdef double cohesion = 0.
        
        for i in range(size-1):
            for j in range(i+1,size):
                cohesion = cohesion + condensed_distance_matrix[cluster[i],cluster[j]]
        
        return weight*cohesion
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

cdef class CythonMirrorCohesionCalculator(object):
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
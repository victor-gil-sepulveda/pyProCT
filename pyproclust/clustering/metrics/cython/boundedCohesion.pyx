# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
'''
Created on 12/07/2012

@author: victor
'''

import numpy
cimport numpy
cimport cython
numpy.import_array()

cdef class CythonBoundedCohesionCalculator(object):
    def __init__(self):
        pass
    
    @cython.boundscheck(False)
    cpdef double evaluate(self,clustering,condensed_distance_matrix)  except *:
        cdef double weight,cohesion,max_cohesion
        cdef int size, i, j, total_number_of_elements
        
        total_number_of_elements = clustering.total_number_of_elements

        if clustering.total_number_of_elements > 0:
            max_cohesion = numpy.sum(condensed_distance_matrix.get_data())/total_number_of_elements
            total_cohesion = 0
            for c in clustering.clusters:
                size = c.get_size()
                weight = 1. / size
                cohesion = 0.
                for i in range(size-1):
                    for j in range(i+1,size):
                        cohesion = cohesion + condensed_distance_matrix[c[i],c[j]]
                total_cohesion +=  weight*cohesion
            return total_cohesion / max_cohesion
        else:
            return 0.
'''
Created on 23/04/2012

@author: victor
'''
import numpy
cimport numpy

import math

DOUBLE = numpy.double
ctypedef numpy.double_t DOUBLE_t

INT = numpy.double
ctypedef numpy.int_t INT_t

def kth_elements_distance(int element, 
                          numpy.ndarray[INT_t] klist, 
                          numpy.ndarray[DOUBLE_t] distances_to_all_other_elements_buffer, 
                          condensed_distance_matrix):
    
    cdef int i
    cdef int k
    
    # Distances of one element vs all the others
    for i in range(condensed_distance_matrix.row_length):
        distances_to_all_other_elements_buffer[i] = condensed_distance_matrix[element,i]
    
    # Pick the kth elements
    distances_to_all_other_elements_buffer.sort(kind='mergesort')
    
    kth_elems = [distances_to_all_other_elements_buffer[k] for k in klist]
    
    return kth_elems

def k_dist(numpy.ndarray[INT_t] klist, 
           numpy.ndarray[DOUBLE_t] buffer, 
           condensed_distance_matrix):
    cdef int N = condensed_distance_matrix.row_length
    cdef int KLEN =  len(klist)
    cdef int i
    
    # for all the elements pick their kth elements
    np_k_dist_matrix = numpy.array([kth_elements_distance(i, klist, buffer,condensed_distance_matrix) for i in range(condensed_distance_matrix.row_length)])

# Now we have:
# element        k_1 k_2 ... K
#    1          [ x   x  ... x ]
#    2          [ x   x  ... x ]
#    ...
#    N          [ x   x  ... x ]
    
# Reshape the matrix
    np_k_dist_matrix = np_k_dist_matrix.T
    
# And now we have
# k         el1  el2 ...  elN
# 20      [  x    x  ...    x ]
# 40      [  x    x  ...    x ]
# ...
# KMAX    [  x    x  ...    x ]

# Rows have to be sorted
    for i in range(KLEN):
        np_k_dist_matrix[i].sort(kind='mergesort')
        
    return np_k_dist_matrix

def dbscan_param_space_search(max_noise, condensed_distance_matrix):
    """
    Does the search of suitable 
    """
    cdef numpy.ndarray[DOUBLE_t] buffer
    cdef int N = condensed_distance_matrix.row_length
    
    # As indicated von Luxburg, 2007) k is in the range log(n)
    klist = k_scale_gen(math.log(N))
    
    buffer = numpy.empty(N)
    kdist_matrix = k_dist(klist,condensed_distance_matrix)

    number_of_elements = condensed_distance_matrix.row_length

    index_of_noise_limit =  int(number_of_elements - (min(max_noise+2.5,100)*0.01*number_of_elements))

    params = []
    
    for i in range(len(klist)):
        params.append((klist[i],kdist_matrix[i][index_of_noise_limit]))
    
    del kdist_matrix
    
    return params
    
def k_scale_gen(max_elements):
    """
    Generates a list of ks as powers of 2.
    """
    k_scale = []
    accum = 1
    try:
        range_max = int(math.ceil(math.log(max_elements,2)))+1
    except ValueError:
        range_max = 1
    
    for i in range(range_max): #@UnusedVariable
        accum *= 2
        k_scale.append(accum)
    return k_scale

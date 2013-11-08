'''
Created on 23/04/2012

@author: victor
'''
import numpy
cimport numpy
import cython
import math

DOUBLE = numpy.double
ctypedef numpy.double_t DOUBLE_t

INT = numpy.double
ctypedef numpy.int_t INT_t

def kth_elements_distance(int element, 
                          numpy.ndarray[INT_t] klist, 
                          numpy.ndarray[DOUBLE_t] distances_to_all_other_elements_buffer, 
                          condensed_distance_matrix):
    
    return cython_kth_elements_distance(element, klist, distances_to_all_other_elements_buffer, condensed_distance_matrix)

cdef numpy.ndarray[DOUBLE_t] cython_kth_elements_distance(int element, 
                          numpy.ndarray[INT_t] klist, 
                          numpy.ndarray[DOUBLE_t] distances_to_all_other_elements_buffer, 
                          condensed_distance_matrix):
    
    cdef int i
    cdef int k
    cdef int KLEN = len(klist)
    cdef int N = condensed_distance_matrix.row_length
    cdef numpy.ndarray[DOUBLE_t] kth_elems = numpy.empty(KLEN)
    
    # Distances of one element vs all the others
    for i in range(N):
        distances_to_all_other_elements_buffer[i] = condensed_distance_matrix[element,i]
    
    # Pick the kth elements
    distances_to_all_other_elements_buffer.sort(kind='quicksort')
    for i in range(KLEN):
        k = klist[i]
        kth_elems[i] = distances_to_all_other_elements_buffer[k]
    
    return kth_elems

@cython.boundscheck(False)
def k_dist(numpy.ndarray[INT_t] klist, 
           numpy.ndarray[DOUBLE_t] buffer, 
           condensed_distance_matrix):
    
    cdef int N = condensed_distance_matrix.row_length
    cdef int KLEN = len(klist)
    cdef int i = 0
    cdef numpy.ndarray[DOUBLE_t, ndim=2] np_k_dist_matrix = numpy.empty((KLEN,N))
    print np_k_dist_matrix
    
    # For all the elements and all Ks pick their kth element distances
    for i in range(N):
        np_k_dist_matrix[:,i] = cython_kth_elements_distance(i, klist, buffer, condensed_distance_matrix) 

    # Now we have
    # k         el1  el2 ...  elN
    # 20      [  x    x  ...    x ]
    # 40      [  x    x  ...    x ]
    # ...
    # KMAX    [  x    x  ...    x ]
    
    # Rows have to be sorted
    for i in range(KLEN):
        np_k_dist_matrix[i].sort(kind='quicksort')
        
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

'''
Created on 23/04/2012

@author: victor
'''
# gcc -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I/usr/lib/python2.7/dist-packages/numpy/core/include -I/usr/include/python2.7 -c cythonDbscanTools.c -o build/temp.linux-x86_64-2.7/cythonDbscanTools.o -O3
# gcc -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions -Wl,-z,relro build/temp.linux-x86_64-2.7/cythonDbscanTools.o -o /home/victor/workspaces/Python/pyProClust/pyproclust/algorithms/dbscan/cython/cythonDbscanTools.so

import numpy
cimport numpy
import cython
import math
from multiprocessing import Pool

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

def dbscan_param_space_search(max_noise, max_eps_tries, number_of_elements, klist, kdist_matrix):
    """
    Does the search of suitable parameters for DBSCAN.
    First it generates a grid of minpts-eps values based on the noise limit imposed by the user.
    Then a new value is added based on (Zhou et al. 2012)
    """

    #MIN_NOISE = 5%
    index_for_min_noise = max(0, int(number_of_elements - 0.05*number_of_elements))
    index_for_max_noise =  int(number_of_elements - (min(max_noise+2.5,100)*0.01*number_of_elements))
    noise_stride = (index_for_min_noise - index_for_max_noise) / max_eps_tries
    
    params = []
    for i in range(index_for_max_noise, index_for_min_noise, noise_stride):
        for j in range(len(klist)):
            params.append((klist[j],kdist_matrix[j][i]))
    
    del kdist_matrix
    
    return params
    
def k_scale_gen(max_elements):
    """
    Generates a list of ks as powers of 2.
    """
    # As indicated von Luxburg, 2007) k is in the range log(n)
    try:
        range_max = int(math.ceil(math.log(max_elements,2)))+1
    except ValueError:
        return []
    
    return numpy.array([2**i for i in range(1,range_max)])

def zhou_adaptative_determination(matrix):
    """
    From Zhou et al. 2012 at Journal of Information and Computational Science
    """
    Eps = matrix.calculateMean()
    Minpts = math.floor(numpy.sum([0]+[len(matrix.element_neighbors_within_radius(i,Eps)) for i in range(matrix.row_length)]) / matrix.row_length)
    return [(Minpts,Eps)]

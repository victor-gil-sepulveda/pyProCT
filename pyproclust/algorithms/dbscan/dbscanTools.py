'''
Created on 23/04/2012

@author: victor
'''
import numpy
import math
import time
import pickle

def kth_elements(element,klist,condensed_distance_matrix):
# Old version
#    all_elements = []
#    for i in range(condensed_distance_matrix.row_length):
#        all_elements.append(condensed_distance_matrix[element,i])
    # Distances of one element vs all the others
    all_elements = numpy.array([condensed_distance_matrix[element,i] for i in range(condensed_distance_matrix.row_length)])
    # Pick the kth elements
    all_elements.sort(kind='mergesort')
    
    kth_elems = []
    for k in klist:
        kth_elems.append(all_elements[k])
    return numpy.array(kth_elems)

def kdist(klist,condensed_distance_matrix):
    k_dist_matrix = []
    
    # for all the elements pick their kth elements
#    t0 = time.time()
    
    for i in range(condensed_distance_matrix.row_length):
        k_dist_matrix.append(kth_elements(i, klist, condensed_distance_matrix))
    
#    print "It took",time.time() - t0, "seconds to calculate the kdists."
# Now we have:
# element        k_1 k_2 ... K
#    1          [ x   x  ... x ]
#    2          [ x   x  ... x ]
#    ...
#    N          [ x   x  ... x ]
    
    # Now reshape the matrix
    np_k_dist_matrix = numpy.array(k_dist_matrix)
    shape = np_k_dist_matrix.shape
    result =  np_k_dist_matrix.reshape((shape[1],shape[0]))
# And now we have
# k         el1  el2 ...  elN
# 20      [  x    x  ...    x ]
# 40      [  x    x  ...    x ]
# ...
# KMAX    [  x    x  ...    x ]

    # But we need the sorted rows!
    for k_array in result:
        k_array.sort(kind='mergesort')
        
    return result


def dbscan_param_space_search(max_noise, condensed_distance_matrix):
    """
    Does the search of suitable 
    """
    num_elements = condensed_distance_matrix.row_length
    
    # As indicated von Luxburg, 2007) k is in the range log(n)
    klist = k_scale_gen(math.log(num_elements))
    
    kdist_matrix = kdist(klist,condensed_distance_matrix)

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
        range_max = int(math.ceil(math.log(max_elements,2)))-1
    except ValueError:
        range_max = 1
    
    for i in range(range_max): #@UnusedVariable
        accum *= 2
        k_scale.append(accum)
    return k_scale

"""
Created on 23/04/2012

@author: victor
"""
import numpy
import math

def kth_elements_distance(element, klist, distances_to_all_other_elements_buffer, condensed_distance_matrix):
    # Distances of one element vs all the others
    for i in range(condensed_distance_matrix.row_length):
        distances_to_all_other_elements_buffer[i] = condensed_distance_matrix[element,i]

    # Pick the kth elements
    distances_to_all_other_elements_buffer.sort(kind='mergesort')

    kth_elems = [distances_to_all_other_elements_buffer[k] for k in klist]

    return kth_elems

def k_dist(klist, buffer, condensed_distance_matrix):

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
    for i in range(len(np_k_dist_matrix)):
        np_k_dist_matrix[i].sort(kind='mergesort')

    return np_k_dist_matrix

def dbscan_param_space_search(max_noise, max_eps_tries, number_of_elements, klist, kdist_matrix):
    """
    Does the search of suitable parameters for DBSCAN.
    First it generates a grid of minpts-eps values based on the noise limit imposed by the user.
    Then a new value is added based on (Zhou et al. 2012)
    """

    #MIN_NOISE = 5%
    index_for_min_noise = max(0, int(number_of_elements - 0.05*number_of_elements))
    index_for_max_noise =  int(number_of_elements - (min(max_noise,100)*0.01*number_of_elements)-1)
    noise_stride = max(1,(index_for_min_noise - index_for_max_noise) / max_eps_tries)

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

def zhou_adaptative_determination(kdist_matrix, matrix):
    """
    From Zhou et al. 2012 at Journal of Information and Computational Science
    """
    N = matrix.row_length
    parameters = []
    Eps_estimations = numpy.mean(kdist_matrix, 1)

    for Eps in Eps_estimations:
        Minpts = math.floor(numpy.sum([0]+[len(matrix.element_neighbors_within_radius(i,Eps)) for i in range(N)]) / N)
        parameters.append((Minpts,Eps))

    return parameters

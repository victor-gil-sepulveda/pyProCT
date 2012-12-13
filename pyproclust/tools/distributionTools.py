'''
Created on 16/05/2012

@author: victor
'''
import numpy

def get_distances_for_elems(elems,center,condensed_distance_matrix):
    """
    
    """
    distances = []
    for e in elems:
        distances.append(condensed_distance_matrix[center,e])
    return distances

def get_distance_std_dev_for_elems(elems,center,condensed_distance_matrix):
    """
    
    """
    distances = []
    for e in elems:
        distances.append(condensed_distance_matrix[center,e])
    return numpy.std(distances)


"""
Created on 16/05/2012

@author: victor
"""
import numpy

def get_distances_for_elems(elems,center,condensed_distance_matrix):
    """
    Extracts the distances of a list of elements to another element using the information of a distance matrix.
    
    @param elems: The array containing the elements we want to calculate the distance.
    @param center: The element to which get the distance. 
    @param condensed_distance_matrix: The distance matrix containing the distance data between elements.
    
    @return: An aray with all the distances.
    """
    distances = []
    for e in elems:
        distances.append(condensed_distance_matrix[center,e])
    return distances

def get_distance_std_dev_for_elems(elems, center, condensed_distance_matrix):
    """
    A wrapper over get_distances_for_elems to get the standard deviation.
    @see: get_distances_for_elems
    
    @param elems: The array containing the elements we want to calculate the distance.
    @param center: The element to which get the distance. 
    @param condensed_distance_matrix: The distance matrix containing the distance data between elements.
    
    @return: The std deviation of the distances.
    """
    return numpy.std(get_distances_for_elems(elems,center,condensed_distance_matrix))


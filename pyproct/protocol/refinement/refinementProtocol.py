'''
Created on 17/09/2012

@author: victor
'''
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.clustering.cluster import Cluster

def pick_all_elements_from_clusters(clusters):
    """
    Puts all the elements of a list of clusters into a list.
    """
    all_elements = []
    for c in clusters:
        all_elements.extend(c.all_elements)
    return all_elements

def redefine_clusters_with_map(clusters,elements_map):
    """
    It renames the elements of a list of clusters using a map array.
    """
    new_clusters = []
    for cluster in clusters:
        new_cluster_elems = []
        for element in cluster.all_elements:
            new_cluster_elems.append(elements_map[element])

        new_clusters.append(Cluster(None, new_cluster_elems))
    return new_clusters

def recreate_matrix(matrix, number_of_elements, all_elements_map):
    """
    Picks a portion of a distance matrix and returns it as another matrix.
    """
    c_matrix = CondensedMatrix([0.0]*int(number_of_elements*(number_of_elements-1)/2.))
    for i in range(number_of_elements):
        ele_i = all_elements_map[i]
        for j in range(number_of_elements):
            ele_j = all_elements_map[j]
            c_matrix[i,j] = matrix[ele_i,ele_j]
    return c_matrix
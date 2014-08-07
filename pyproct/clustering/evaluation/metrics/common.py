"""
Created on 06/06/2013

@author: victor
"""
def get_intra_cluster_distances(cluster, matrix):
    distances = []
    cluster_elements = cluster.all_elements
    for i in range(len(cluster.all_elements)-1):
        for j in range(i+1,len(cluster.all_elements)):
            distances.append(matrix[cluster_elements[i],cluster_elements[j]])
    return distances

def get_inter_cluster_distances(i, j, clusters, matrix):
    distances = []
    for cluster_i_element in clusters[i].all_elements:
        for cluster_j_element in clusters[j].all_elements:
            distances.append(matrix[cluster_i_element, cluster_j_element])
    return distances

def get_inter_cluster_prototype_distances(clusters, matrix):
    """
    Precondition: cluster medoids have been updated.
    """
    distances = []
    for i in range(len(clusters)-1):
        for j in range(i+1,len(clusters)):
            distances.append(matrix[clusters[i].prototype,clusters[j].prototype]) 
    return distances

def get_distances_of_elements_to(from_this, to_those, matrix):
    distances = []
    for element in to_those:
        distances.append(matrix[from_this,element])
    return distances

def update_medoids(clustering, matrix):
    for cluster in clustering.clusters:
        cluster.prototype = cluster.calculate_medoid(matrix)
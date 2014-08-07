'''
Created on 27/02/2014

@author: victor
'''
import numpy
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.evaluation.metrics.common import get_distances_of_elements_to

def calculate_mean_center_differences(decomposed_cluster, matrix):
    """
    Given a mixed decomposed cluster, it calculates the mean of all center differences (giving a qualitative
    view of how separated the inner subclusters are).

    @param decomposed_cluster: A MIXED decomposed cluster.
    @param matrix: The condensed distance matrix used.

    @return: The mean of center distances.
    """
    centers = []
    for traj_id in decomposed_cluster:
        cluster = Cluster(None, decomposed_cluster[traj_id])
        centers.append(cluster.calculate_medoid(matrix))
    center_distances = []
    for i in range(len(centers)-1):
        for j in range(i+1, len(centers)):
            center_distances.append(matrix[centers[i],centers[j]])
    return numpy.mean(center_distances)


def calculate_distance_stats(elements, matrix):
    """
    Calculates the mean, dispersion and radius of all the distances to the central element of a set of
    elements.

    @param elements: The elements we are working with.
    @param matrix: The used condensed matrix.

    @return: Mean, std deviation and radius of all the elements with respect to their central element.
    """
    cluster = Cluster(None, elements)
    medoid = cluster.calculate_medoid(matrix)

    # We also get a 0 distance from the medoid vs itself (it is contained in 'elements')
    distances = get_distances_of_elements_to(medoid, elements, matrix)
    return numpy.mean(distances), numpy.std(distances), numpy.max(distances)


def getAllElements(decomposed_cluster):
    """
    Returns a list of all the elements of a decomposed cluster (the elements of the original cluster)

    @param decomposed_cluster: One decomposed cluster (trajectory_id::elements)

    @return: An array containing all the elements of the decomposed cluster.
    """
    all_elements = []
    for traj_id in decomposed_cluster:
        all_elements.extend(decomposed_cluster[traj_id])
    return all_elements


def mergeSeparatedClusters(separated_decomposed_clusters):
    """
    Returns a list of all the clusters of a separated decomposed clustering.

    @param separated_decomposed_clusters: One separated decomposed clustering (pure/mixed::cluster_id::traj_id::elements)

    @return: An array containing all the decomposed clusters of the separated clusterings.
    """
    decomposed_clusters = []
    for ctype in separated_decomposed_clusters:
        for cluster_id in separated_decomposed_clusters[ctype]:
            decomposed_clusters.append(separated_decomposed_clusters[ctype][cluster_id])
    return decomposed_clusters
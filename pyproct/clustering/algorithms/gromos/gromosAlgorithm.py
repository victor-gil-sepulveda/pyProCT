"""
Created on 14/02/2012

@author: victor
"""
from pyproct.clustering.cluster import cluster_from_tuple
from pyproct.clustering.algorithms.gromos.gromosAlgorithmTools import eliminate_cluster_from_node_list
from pyproct.clustering.clustering import Clustering
import sys
import numpy


class GromosAlgorithm(object):
    """
    TODOC
    """
    def __init__(self,condensed_matrix, **kwargs):
        self.condensed_matrix = condensed_matrix

    def perform_clustering(self,kwargs):
        """
        Does the actual clustering.
        """
        cutoff = kwargs["cutoff"]

        try:
            max_clusters = kwargs["max_clusters"]
        except KeyError:
            max_clusters = sys.maxint

        nodes = range(self.condensed_matrix.row_length)
        clusters = []
        elements_already_clustered = 0
        iteration = 0
        # Do it while there are nodes left
        while not len(nodes) == 0 and not len(clusters) >= max_clusters:
            cluster = self.__do_one_iteration(nodes,cutoff)
            clusters.append(cluster)
            elements_already_clustered = elements_already_clustered + cluster.get_size()
            if elements_already_clustered + len(nodes) > self.condensed_matrix.row_length:
                print "[ERROR :: GROMOS perform_clustering] ", elements_already_clustered + len(nodes), iteration
                exit(1)
            iteration = iteration + 1
            
        return Clustering(clusters,details="GROMOS (cutoff = "+str(cutoff)+")")

    def __do_one_iteration(self,remaining_nodes,cutoff):
        """
        IMPORTANT: Prototype and medoid may not be the same element (ex @ http://www.jstatsoft.org/v01/i04/paper).
        """
        
        (node,len_neigh) = self.condensed_matrix.choose_node_with_higher_cardinality(remaining_nodes,cutoff) #@UnusedVariable
        neighbours = self.condensed_matrix.get_neighbors_for_node(node,remaining_nodes,cutoff)
        cluster_tuple =  (node, neighbours)
        cluster = cluster_from_tuple(cluster_tuple)
        eliminate_cluster_from_node_list(remaining_nodes,cluster)
        return cluster




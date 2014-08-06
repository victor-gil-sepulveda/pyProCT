"""
Created on 19/04/2012

@author: victor
"""
import random
from pyproct.clustering.cluster import gen_clusters_from_class_list
from pyproct.clustering.clustering import Clustering

class RandomClusteringAlgorithm(object):

    def __init__(self,condensed_matrix, **kwargs):
        """
        Constructor.
        """
        self.condensed_matrix = condensed_matrix

    def perform_clustering(self, kwargs):
        """
        Creates a clustering where the clusters have been created by random selection of
        the elements in the dataset. It will create a random number of clusters if "max_num_of_clusters" is
        defined, or an exact number of clusters this clusters if "num_clusters" is defined.
        """
        num_of_nodes = self.condensed_matrix.row_length
        num_of_clusters = 0
        elements_per_cluster = 0
        max_num_of_clusters = 0
        try:
            num_of_clusters = kwargs["num_clusters"]
            elements_per_cluster = max(1, num_of_nodes / num_of_clusters)
        except KeyError:
            try:
                max_num_of_clusters = kwargs["max_num_of_clusters"]
                num_of_clusters = random.randint(1, max_num_of_clusters)
            except KeyError:
                num_of_clusters = 2
                elements_per_cluster = 1

        node_class = []

        try:
            elements_per_cluster = max(1, num_of_nodes / num_of_clusters)
        except:
            elements_per_cluster = 1

        for i in range(num_of_clusters):
            node_class.extend([i]*elements_per_cluster)

        while len(node_class) < num_of_nodes:
            node_class.append(0)

        random.seed()
        random.shuffle(node_class)

        clusters = gen_clusters_from_class_list(node_class)

        return Clustering(clusters, details = "Random (max_num_of_clusters = %d, num_of_clusters = %d )"%(max_num_of_clusters, num_of_clusters))
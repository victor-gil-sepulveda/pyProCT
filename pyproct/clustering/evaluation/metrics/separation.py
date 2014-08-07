"""
Created on 06/06/2013

@author: victor
"""
"""
Created on 09/01/2013

@author: victor
"""
from pyproct.clustering.evaluation.metrics.cohesion import CohesionCalculator

class SeparationCalculator(object):

    def __init__(self):
        pass

    def evaluate(self, clustering, condensed_distance_matrix, cohesions = None):
        all_clusters = clustering.clusters

        my_cohesions_dic = {}
        if (cohesions is None):
            #Calculate_cohesions if not given
            cohesion_calctor = CohesionCalculator()
            for cluster in all_clusters:
                my_cohesions_dic[cluster] =  cohesion_calctor.evaluate_cluster(cluster, condensed_distance_matrix)
        else:
            for i, cluster in enumerate(all_clusters):
                my_cohesions_dic[cluster] = cohesions[i]

        total_separation = 0
        for i, cluster in enumerate(all_clusters):
            separation = self.cluster_separation(cluster, clustering, my_cohesions_dic[cluster], condensed_distance_matrix)
            total_separation = total_separation + separation

        return total_separation

    def cluster_separation(self, cluster, clustering, clustering_cohesion, condensed_distance_matrix):
        """
        Returns the cohesion plus separation value of a cluster. The condensed matrix
        given as parameter stores the distances of the elements of the dataset used to
        extract the cluster.
        """
        if clustering_cohesion > 0:
            weight = 1./clustering_cohesion
            sep_and_cohe = 0.0
            ## I'm inside?
            where_am_i = clustering.cluster_index(cluster)

            for i in range(len(clustering.clusters)):
                if i != where_am_i :
                    c_j = clustering.clusters[i]
                    sep_and_cohe = sep_and_cohe + self.__between_cluster_distance(cluster,c_j,condensed_distance_matrix)
            return weight*sep_and_cohe
        else:
            return 0. # not defined in this case, could be numpy.finfo(numpy.float32).max

    def __between_cluster_distance(self,cluster_1,cluster_2,condensed_distance_matrix):
        """
        Calculates the 'cohesion' of one cluster vs other.
        Precondition: Clusters don't have shared elements.
        """
        mixed_cohesion = 0
        for c_i in cluster_1.all_elements:
            for c_j in cluster_2.all_elements:
                mixed_cohesion = mixed_cohesion + condensed_distance_matrix[c_i,c_j]
        return mixed_cohesion

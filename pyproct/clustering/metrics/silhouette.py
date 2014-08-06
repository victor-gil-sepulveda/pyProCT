"""
Created on 02/05/2012

@author: victor
"""
import numpy

class SilhouetteCoefficientCalculator(object):
    """
    Between -1 and 1. 1 Is the value for the better clustering (best separated and with best cohesion.
    0 can be a not-so-good-and-not-so-bad clustering, but a random clustering would get results
    around this.
    """
    def __init__(self):
        self.warn_misclustered_elements = False

    def evaluate(self, clustering, condensed_distance_matrix, element = None, cluster = None):
        """
        Evaluates the silhouette factor for 3 types of elements in the clusterization
        hierarchy:
        - element silhouette: it needs to define the element and the cluster to which
        this element belongs
        - cluster: only needs to define the cluster.
        - clusterization: calculates the silhouette factor of a whole clusterization.
        """
        self.warn_misclustered_elements = False

        if not element and not cluster and not clustering:
            print "[Error SilhouetteCoefficientCalculator:evaluate] you may use at least one parameter"
            exit(-1)

        if (element and clustering) or\
           (cluster and clustering) or\
           (element and not cluster) or\
           (element and cluster and clustering):
            print "[Error SilhouetteCoefficientCalculator:evaluate] wrong parametrization"
            exit(-1)

        if element:
            result = self.__one_element_silhouette(element,cluster,condensed_distance_matrix)

        if cluster:
            result = numpy.mean(self.__one_cluster_partial_silhouette(cluster,clustering,condensed_distance_matrix))

        if clustering:
            result = numpy.mean(self.__one_clusterization_partial_silhouette(clustering,condensed_distance_matrix))

        if self.warn_misclustered_elements:
            print "[WARNING SilhouetteCoefficientCalculator:evaluate] Some elements may be misclassified and Silhouette was going to return NaN."

        return result

    def __one_clusterization_partial_silhouette(self,clusterization,condensed_distance_matrix):
        """
        Calculates the partial results of the silhouette coefficient for a clusterization.
        """
        cluster_silhouettes = []
        for c in clusterization.clusters:
            cluster_silhouettes.extend(self.__one_cluster_partial_silhouette(c,clusterization,condensed_distance_matrix))
        return cluster_silhouettes

    def __one_cluster_partial_silhouette(self,cluster,clusterization,condensed_distance_matrix):
        """
        Calculates the partial results of the silhouette coefficient for a cluster.
        """
        silhouette_factors=[]
        for element in cluster:
            silhouette_factors.append(self.__one_element_silhouette(element,cluster,clusterization,condensed_distance_matrix))
        return silhouette_factors

    def __one_element_silhouette(self,element,cluster,clusterization,condensed_distance_matrix):
        """
        element is inside cluster
        """
        a_i = self.__get_average_distance_with_my_cluster(element,cluster,condensed_distance_matrix)
        where_am_i = clusterization.cluster_index(cluster)
        b_i = numpy.finfo(numpy.float32).max
        for i in range(len(clusterization.clusters)):
            if where_am_i != i:
                b_i = min(b_i,self.__get_average_distance_with_cluster(element,clusterization.clusters[i],condensed_distance_matrix))
        if b_i == 0:
            self.warn_misclustered_elements = True # if b1 equals 0, means that we did something wrong while clustering (an element has found a cluster
            # where all other elements are equal and equal to itself, this element is incorrectly placed)
            return -1
        else:
            return (b_i-a_i)/max(a_i,b_i)

    def __sum_cluster_distances(self,element,cluster,condensed_distance_matrix):
        """
        Returns the sum of the distances of one element vs one cluster.
        """
        distance_sum = 0.0
        for e in cluster.all_elements:
            distance_sum += condensed_distance_matrix[e,element]
        return distance_sum

    def __get_average_distance_with_my_cluster(self,element,cluster,condensed_distance_matrix):
        """
        Calculates the average distance of one element to all the elements of his cluster.
        """
        if cluster.get_size() == 1:
            # Then we will only rely in b_i, the bigger is b_i, the bigger the coefficient (it
            # has sense, as b_i is calculating the separation, and the separation is better
            # as it becomes bigger
            return 1 ## a_i is non negative and 0 is the number leading to the bigger coef.
        else:
            distance_sum = self.__sum_cluster_distances(element,cluster,condensed_distance_matrix)
            return distance_sum / (cluster.get_size()-1)

    def __get_average_distance_with_cluster(self,element,cluster,condensed_distance_matrix):
        """
        Calculates the average distance of one element to all the elements of one cluster.
        """
        distance_sum = self.__sum_cluster_distances(element,cluster,condensed_distance_matrix)
        return distance_sum / cluster.get_size()

'''
Created on 06/06/2013

@author: victor
'''
import numpy

class DunnCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        dmin = numpy.min([self.min_intracluster_distance(c) for c in clustering.clusters])
        dmax = self.max_intercluster_distance(clustering, matrix)
        dunn_index = dmin / dmax
        return dunn_index
        
    @classmethod
    def min_intracluster_distance(cls, cluster, matrix):
        distances = []
        for i in range(len(cluster.all_elements)-1):
            for j in range(i+1,len(cluster.all_elements)):
                distances.append(matrix[i,j])
        return numpy.min(distances)
    
    @classmethod
    def max_intercluster_distance(cls, clustering, matrix):
        distances = []
        for i in range(len(clustering.clusters)-1):
            for j in range(i+1,len(clustering.clusters)):
                distances.append(cls.max_pairwise_intercluster_distance(i, j, clustering.clusters, matrix))
        return numpy.max(distances)
    
    @classmethod
    def max_pairwise_intercluster_distance(cls, i, j, clusters, matrix):
        distances = []
        for cluster_i_element in clusters[i].all_elements:
            for cluster_j_element in clusters[j].all_elements:
                distances.append(matrix[cluster_i_element, cluster_j_element])
        return numpy.max(distances)
    
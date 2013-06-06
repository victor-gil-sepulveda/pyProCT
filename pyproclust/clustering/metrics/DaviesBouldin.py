'''
Created on May 30, 2013

@author: victor
'''
import numpy

class DaviesBouldinCalculator(object):
    def __init__(self):
        pass
    
    def evaluate(self, clustering, matrix):
        
        # Medoids will be the new prototypes
        self.update_medoids(clustering, matrix)
        
        # We calculate the average distances to the medoid for each of the clusters
        average_distances = [self.calculate_average_distance_from_prototype(cluster, matrix) for cluster in clustering.clusters]
        
        # Then we can do the actual calculation
        db_index = self.db_index_calculation( average_distances, clustering.clusters, matrix)
        
        return db_index
    
    @classmethod
    def update_medoids(cls, clustering, matrix):
        for cluster in clustering.clusters:
            cluster.prototype = cluster.calculate_medoid(matrix)
    
    @classmethod
    def db_index_calculation(cls, d, clusters, matrix):
        db_index = 0
        for i in range(len(clusters)-1):
            db_index += cls.max_db_term(i, d, clusters, matrix)
        return db_index/len(clusters)
        
    @classmethod
    def max_db_term(cls, i, d, clusters, matrix):
        terms = []
        for j in range(i+1, len(clusters)):
            if i != j:
                terms.append(cls.d_b_term_calculation(i,j,d,clusters, matrix))
    
    @classmethod
    def d_b_term_calculation(cls, i, j, d, clusters, matrix):
        c_i, c_j = clusters[i].prototype, clusters[j].prototype
        return (d[i]+d[j]) / (matrix[c_i, c_j])
    
    @classmethod    
    def calculate_average_distance_from_prototype(cls, cluster, matrix):
        distances = []
        proto = cluster.prototype
        for element in cluster.all_elements:
            distances.append(matrix[element][proto])
        return numpy.mean(distances)   
    
        
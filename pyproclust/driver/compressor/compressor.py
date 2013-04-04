'''
Created on 04/04/2013

@author: victor
'''
import pyproclust.protocol.saveTools as saveTools
import math
from pyproclust.algorithms.kmedoids import kMedoidsAlgorithm
from pyproclust.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm

class Compressor(object):

    def __init__(self, parameters):
        self.type = parameters["type"]
        self.parameters = parameters
        
    def compress(self, clustering, pdb_name, workspace_handler, trajectory_handler, matrix_handler):
        representatives = []
        if self.compression_type == "RANDOM":
            representatives = self.__naive_compression(clustering, matrix_handler)
            
        elif self.compression_type == "KMEANS":
            representatives = self.__kmeans_compression()
        
        else:
            print "[ERROR Compressor::compress] The compression type does not exist (%s)"%(self.type)
        saveTools.save_representatives(representatives, pdb_name, workspace_handler, trajectory_handler)
    
    def __naive_compression(self, clustering, matrix_handler):
        """
        
        @return: The nth most representative elements of a clustering.
        """
        number_of_final_frames = self.parameters["final_number_of_frames"]
        return  clustering.get_proportional_size_representatives(number_of_final_frames, 
                                                                matrix_handler.distance_matrix)
        
    def __kmeans_compression(self, clustering, matrix_handler):
        """
        """
        representatives = []
        for cluster in clustering.clusters:
            cluster_size = float(cluster.get_size())
            expected_cluster_elements = cluster_size / (clustering.total_number_of_elements * self.parameters["final_number_of_frames"])
            expected_cluster_elements = int(math.ceil(expected_cluster_elements))
            kmedoids = KMedoidsAlgorithm(matrix_handler.distance_matrix)
            clustering = kmedoids.perform_clustering({
                                                      "k":expected_cluster_elements,
                                                      "seeding_type":"RANDOM"
                                                      })
            representatives.extend(clustering.get_medoids(matrix_handler.distance_matrix))
            
        
            
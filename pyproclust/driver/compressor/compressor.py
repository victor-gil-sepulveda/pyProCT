'''
Created on 04/04/2013

@author: victor
'''
import pyproclust.protocol.saveTools as saveTools
import math
from pyproclust.algorithms.kmedoids import kMedoidsAlgorithm
from pyproclust.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproclust.protocol.refinement.refinementProtocol import redefine_clusters_with_map,\
    recreate_matrix
from pyproclust.clustering.cluster import Cluster

class Compressor(object):

    def __init__(self, parameters):
        
        self.parameters = parameters
        
    def compress(self, clustering, pdb_name, workspace_handler, trajectory_handler, matrix_handler):
        representatives = []
        compression_type = self.parameters["type"]
        
        if compression_type == "RANDOM":
            representatives = self.__naive_compression(clustering, matrix_handler)
        
        elif compression_type == "KMEDOIDS":
            representatives = self.__kmedoids_compression(clustering, matrix_handler)
        
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
        
    def __kmedoids_compression(self, clustering, matrix_handler):
        """
        """
        representatives = []
        for cluster in clustering.clusters:
            # Guess 'correct' number of elements for this cluster
            cluster_size = cluster.get_size()
            expected_cluster_elements = cluster_size * (float(self.parameters["final_number_of_frames"]) / clustering.total_number_of_elements)
            expected_cluster_elements = int(math.ceil(expected_cluster_elements))
            
            # Remap this cluster elements to get a new matrix
            elements_map = []
            for e in cluster.all_elements:
                elements_map.append(e) # Position [0] for element x1 etc... 
            remapped_matrix = recreate_matrix(matrix_handler.distance_matrix, 
                                              cluster_size, 
                                              elements_map)            
           
            # Prepare and run kmedoids algorithm            
            kmedoids = KMedoidsAlgorithm(remapped_matrix)
#             print "KMEDOIDS:: EXPECTED", expected_cluster_elements, cluster_size, clustering.total_number_of_elements, self.parameters["final_number_of_frames"]
            new_clustering = kmedoids.perform_clustering({
                                                      "k": expected_cluster_elements,
                                                      "seeding_type": "RANDOM"
            })
            
#             print "NEW CLUSTERING SIZE  clusters: %d  elements: %d"%(len(new_clustering.clusters), new_clustering.total_number_of_elements)
            
            remapped_representatives = new_clustering.get_medoids(remapped_matrix)
            fake_cluster = Cluster(None, remapped_representatives)
            
            representatives.extend(redefine_clusters_with_map([fake_cluster], elements_map)[0].all_elements)
        
        return representatives
        
            
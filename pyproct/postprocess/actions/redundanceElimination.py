"""
Created on 04/04/2013

@author: victor
"""
import os
import math
from pyproct.clustering.cluster import Cluster
from pyproct.tools.matrixTools import get_submatrix
from pyproct.clustering.protocol.refinement.Refiner import Refiner
from pyproct.clustering.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproct.postprocess.actions.representatives import save_cluster_elements


class RedundanceEliminationPostAction(object):
    KEYWORD = "compression"

    def __init__(self):
        pass

    def run(self, clustering, postprocessing_parameters, data_handler, workspaceHandler, matrixHandler, generatedFiles):
        compressor = RedundanceElimination(postprocessing_parameters)
        compressed_file_path = compressor.compress(clustering,
                                                   workspaceHandler,
                                                   data_handler,
                                                   matrixHandler,
                                                   postprocessing_parameters)
        generatedFiles.append({
                               "description":"Compressed file",
                               "path":os.path.abspath(compressed_file_path),
                               "type":"pdb"
        })


class RedundanceElimination(object):

    def __init__(self, parameters):
        self.parameters = parameters

    def compress(self, clustering, workspace_handler, data_handler, matrix_handler, options):
        representatives = []
        compression_type = self.parameters.get_value("type", default_value = "KMEDOIDS")

        pdb_name =  self.parameters.get_value("filename", default_value = "compressed")
        pdb_path = os.path.join( workspace_handler["results"],"%s.pdb"%pdb_name)

        if compression_type == "RANDOM":
            representatives = self.__naive_compression(clustering, matrix_handler)

        elif compression_type == "KMEDOIDS":
            representatives = self.__kmedoids_compression(clustering, matrix_handler)

        else:
            print "[ERROR Compressor::compress] The compression type does not exist (%s)"%(self.type)
            exit() 
        
        save_cluster_elements(representatives,
                              pdb_path,
                              data_handler,
                              options)
        return pdb_name

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

            remapped_matrix = get_submatrix(matrix_handler.distance_matrix, cluster.all_elements)

            # Prepare and run kmedoids algorithm
            kmedoids = KMedoidsAlgorithm(remapped_matrix)
#             print "KMEDOIDS:: EXPECTED", expected_cluster_elements, cluster_size, clustering.total_number_of_elements, self.parameters["final_number_of_frames"]
            new_clustering = kmedoids.perform_clustering({
                                                      "k": expected_cluster_elements,
                                                      "seeding_type": "EQUIDISTANT"
            })

#             print "NEW CLUSTERING SIZE  clusters: %d  elements: %d"%(len(new_clustering.clusters), new_clustering.total_number_of_elements)

            # reverse the remapping and add it to representatives
            remapped_representatives = new_clustering.get_medoids(remapped_matrix)
            fake_cluster = Cluster(None, remapped_representatives)

            representatives.extend(Refiner.redefine_cluster_with_map(cluster, fake_cluster).all_elements)

        return representatives


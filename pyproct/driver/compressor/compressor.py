'''
Created on 04/04/2013

@author: victor
'''
import pyproct.protocol.saveTools as saveTools
import math
from pyproct.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproct.clustering.cluster import Cluster
from pyproct.protocol.refinement.Refiner import Refiner
from pyproct.tools.matrixTools import get_submatrix

class Compressor(object):

    def __init__(self, parameters):
        self.parameters = parameters

    def compress(self, clustering, workspace_handler, trajectory_handler, matrix_handler):
        representatives = []
        compression_type = self.parameters["type"]

        pdb_name =  self.parameters['file'] if 'file' in self.parameters else "compressed_pdb"

        if compression_type == "RANDOM":
            representatives = self.__naive_compression(clustering, matrix_handler)

        elif compression_type == "KMEDOIDS":
            representatives = self.__kmedoids_compression(clustering, matrix_handler)

        else:
            print "[ERROR Compressor::compress] The compression type does not exist (%s)"%(self.type)

        return saveTools.save_representatives(representatives,
                                              pdb_name,
                                              workspace_handler,
                                              trajectory_handler,
                                              do_merged_files_have_correlative_models=True,
                                              write_frame_number_instead_of_correlative_model_number=False,
                                              keep_remarks = (lambda params: params['keep_remarks'] if 'keep_remarks' in params else False)(self.parameters))

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


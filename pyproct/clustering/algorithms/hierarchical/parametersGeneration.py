"""
Created on 27/05/2013

@author: victor
"""
from pyproct.clustering.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
import pyproct.clustering.algorithms.hierarchical.hierarchicalTools as  hierarchicalTools

class ParametersGenerator(object):
    HIERARCHICAL_REFINEMENT_VALUE = 200

    def __init__(self, parameters, matrix_handler):
        """
        Class creator.

        @param parameters: Script parameters.

        @param distance_matrix: The distance matrix we are using.
        """
        self.distance_matrix = matrix_handler.distance_matrix
        self.parameters = parameters

    @classmethod
    def get_base_parameters(cls):
        """
        Defines the base parameters needed for each of the algorithms. Each parameter created will be based
        on one of those and must not have more keys than these.

        @return: A dictionary with the base parameters for this algorithm.
        """
        return {
                 "cutoff": None,
                 "method": None
        }

    def get_parameters(self):
        """
        This function creates some parameters to be used with the Hierarchical algorithm.
        @return: A tuple with the generated parameters and a list with its corresponding clusterings
        (in this case parameters and clusterings are obtained at the same time).
        """
        run_parameters = []
        max_clusters = self.parameters["clustering"]["evaluation"]["maximum_clusters"]
        min_clusters = self.parameters["clustering"]["evaluation"]["minimum_clusters"]
        hierarchicalAlgorithm = HierarchicalClusteringAlgorithm(self.distance_matrix)
        clusters_and_cutoff = hierarchicalTools.get_clusters_with_ranged_search(
                                                                        hierarchicalAlgorithm,
                                                                        0.,
                                                                        self.distance_matrix.calculateMean(),
                                                                        min_clusters,
                                                                        max_clusters,
                                                                        ParametersGenerator.HIERARCHICAL_REFINEMENT_VALUE)
        clusterings = []
        cutoffs = []
        for numclusters in clusters_and_cutoff:
            clusterings.append(clusters_and_cutoff[numclusters][1])
            cutoffs.append(clusters_and_cutoff[numclusters][0])

        for cutoff in cutoffs:
            run_parameter = ParametersGenerator.get_base_parameters()
            run_parameter["method"] = 'complete'
            run_parameter["cutoff"] = cutoff
            run_parameters.append(run_parameter)

        return run_parameters, clusterings

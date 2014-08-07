"""
Created on 27/05/2013

@author: victor
"""
import pyproct.clustering.algorithms.dbscan.cython.cythonDbscanTools as dbscanTools
from pyproct.clustering.algorithms.dbscan.cython.cythonDbscanTools import k_scale_gen,\
    k_dist, zhou_adaptative_determination
import numpy
import math

class ParametersGenerator(object):

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
        Defines the base parameters needed for each of the algorithms. Each created parameter will be based
        on one of those and must not have more keys than these.

        @return: A dictionary with the base parameters for this algorithm.
        """
        return {
                 "eps": None,
                 "minpts": None
        }

    def get_parameters(self):
        """
        This function creates some parameters to be used with DBScan.
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []

        # (minpts, eps tuples)
        if "max" in self.parameters["clustering"]["algorithms"]["dbscan"]:
            max_eps_tries = self.parameters["clustering"]["algorithms"]["dbscan"]["max"]
        else:
            max_eps_tries = 10

        num_elements = self.distance_matrix.row_length
        klist = k_scale_gen(math.log(num_elements))
        buffer = numpy.empty(num_elements)
        kdist_matrix = k_dist(klist, buffer, self.distance_matrix)

        dbscan_param_pairs = dbscanTools.dbscan_param_space_search(self.parameters["clustering"]["evaluation"]["maximum_noise"],
                                                                   max_eps_tries,
                                                                   num_elements,
                                                                   klist,
                                                                   kdist_matrix) +\
                            zhou_adaptative_determination(kdist_matrix, self.distance_matrix)

        for (minpts, eps) in dbscan_param_pairs:
            run_parameter = ParametersGenerator.get_base_parameters()
            run_parameter["minpts"] = minpts
            run_parameter["eps"] = eps
            run_parameters.append(run_parameter)

        return run_parameters, []


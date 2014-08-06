"""
Created on 27/05/2013

@author: victor
"""
import math
import hashlib


class ParametersGenerator(object):

    MAX_CLUSTERING_PARAMETERS = 15.

    def __init__(self, parameters, matrix_handler):
        """
        Class creator.

        @param parameters: Script parameters.

        @param distance_matrix: The distance matrix we are using.
        """
        self.parameters = parameters

        if "max" in parameters["clustering"]["algorithms"]["kmedoids"]:
            max_gen_clusterings = float(parameters["clustering"]["algorithms"]["kmedoids"]["max"])
        else:
            max_gen_clusterings = ParametersGenerator.MAX_CLUSTERING_PARAMETERS

        self.num_clusters_step = int(math.ceil((parameters["clustering"]["evaluation"]["maximum_clusters"] - parameters["clustering"]["evaluation"]["minimum_clusters"]) / max_gen_clusterings))
        if self.num_clusters_step < 1:
            self.num_clusters_step = 1

    @classmethod
    def get_base_parameters(cls):
        """
        Defines the base parameters needed for each of the algorithms. Each parameter created will be based
        on one of those and must not have more keys than these.

        @return: A dictionary with the base parameters for this algorithm.
        """
        return {
                 "k": None,
                 "seeding_type": None,
                 "seeding_max_cutoff": None
        }

    def get_parameters(self):
        """
        This function creates some parameters to be used with K-Medoids.
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []
        max_clusters = self.parameters["clustering"]["evaluation"]["maximum_clusters"]
        min_clusters = self.parameters["clustering"]["evaluation"]["minimum_clusters"]
        sizes = range(min_clusters,max_clusters+1,self.num_clusters_step)

        # Defaults
        if not "seeding_type" in self.parameters["clustering"]["algorithms"]["kmedoids"]:
            self.parameters["clustering"]["algorithms"]["kmedoids"]["seeding_type"] = "EQUIDISTANT"

        # This could be subclassed
        if self.parameters["clustering"]["algorithms"]["kmedoids"]["seeding_type"] == "RANDOM":
            tries = self.parameters["clustering"]["algorithms"]["kmedoids"]["tries"] if "tries" in self.parameters["clustering"]["algorithms"]["kmedoids"] else 10

            for i in range(tries):
                for one_size in sizes:
                    run_parameter = ParametersGenerator.get_base_parameters()
                    run_parameter["k"]  = one_size
                    run_parameter["seeding_type"] = "RANDOM"
                    run_parameter["rand_seed"] = int(hashlib.sha256(str({"dic":self.parameters,"iter":i})).hexdigest(),base=16)
                    run_parameters.append(run_parameter)
        else:
            for one_size in sizes:
                run_parameter = ParametersGenerator.get_base_parameters()
                run_parameter["k"]  = one_size
                run_parameter["seeding_type"] = "EQUIDISTANT"
                run_parameters.append(run_parameter)

        return run_parameters, []

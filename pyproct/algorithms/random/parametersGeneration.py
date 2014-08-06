"""
Created on 27/05/2013

@author: victor
"""
class ParametersGenerator(object):

    CLUSTERING_SIZE_STEP = 2

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
                  "num_clusters": None
        }

    def get_parameters(self):
        """
        This function creates some parameters to be used with the random algorithm.
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []
        max_clusters = self.parameters["clustering"]["evaluation"]["maximum_clusters"]
        min_clusters = self.parameters["clustering"]["evaluation"]["minimum_clusters"]
        sizes = range(min_clusters,max_clusters+1,ParametersGenerator.CLUSTERING_SIZE_STEP)
        for one_size in sizes:
            run_parameter = ParametersGenerator.get_base_parameters()
            run_parameter["num_clusters"]  = one_size
            run_parameters.append(run_parameter)

        return run_parameters, []

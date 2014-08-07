"""
Created on 27/05/2013

@author: victor
"""
import numpy

class ParametersGenerator(object):
    GROMOS_DEFAULT_MAX = 25

    def __init__(self, parameters, matrix_handler):
        """
        Class creator.

        @param parameters: Script parameters.

        @param distance_matrix: The distance matrix we are using.
        """
        self.distance_matrix = matrix_handler.distance_matrix
        self.parameters = parameters
        if "max" in parameters["clustering"]["algorithms"]["gromos"]:
            self.max_gen_clusterings = float(parameters["clustering"]["algorithms"]["gromos"]["max"])
        else:
            self.max_gen_clusterings = ParametersGenerator.GROMOS_DEFAULT_MAX
    @classmethod
    def get_base_parameters(cls):
        """
        Defines the base parameters needed for each of the algorithms. Each parameter created will be based
        on one of those and must not have more keys than these.

        @return: A dictionary with the base parameters for this algorithm.
        """
        return {
                  "cutoff": None
        }

    def get_most_separated_elements(self):
        """
        Chooses the 2 most separated elements of the dataset based on the distance stored in distance matrix.
        @return: A dictionary with the elements and their distance value.
        """
        max_s = {
                 "elements":(0,0),
                 "value":0
                 }

        for i in range(self.distance_matrix.row_length):
            for j in range(i+1,self.distance_matrix.row_length):
                if self.distance_matrix[i,j] > max_s["value"]:
                    max_s["elements"] = (i,j)
                    max_s["value"] = self.distance_matrix[i,j]

        return max_s

    def get_most_separated_from_two_elements(self, first, second):
        """
        Gets the element which distance to two elements is maximum.

        @param first: First element to check against.
        @param second: Second element to check against.

        @return: A dictionary with the number of element and the
        """
        max_s = {
                 "element":0,
                 "mean_value":0,
                 "value":0
                 }

        for i in range(self.distance_matrix.row_length):
            my_mean = numpy.mean([self.distance_matrix[i,first],self.distance_matrix[i,second]])
            if my_mean > max_s["mean_value"]:
                max_s["element"] = i
                max_s["mean_value"] = my_mean
                max_s["value"] = numpy.max([self.distance_matrix[i,first],self.distance_matrix[i,second]])
        return max_s

    def get_parameters(self):
        """
        This function creates some parameters to be used with Gromos.
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        most_separated = self.get_most_separated_elements()
        central_element =  self.get_most_separated_from_two_elements(most_separated["elements"][0],most_separated["elements"][1])

        step = central_element["value"] / self.max_gen_clusterings
        run_parameters = []
        cutoff = step
        while cutoff < central_element["value"]:
            run_parameter = ParametersGenerator.get_base_parameters()
            run_parameter["cutoff"] = cutoff
            cutoff += step
            run_parameters.append(run_parameter)

        return run_parameters, []


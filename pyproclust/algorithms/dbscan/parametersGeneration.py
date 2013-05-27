'''
Created on 27/05/2013

@author: victor
'''
import pyproclust.algorithms.dbscan.dbscanTools as dbscanTools

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
        Defines the base parameters needed for each of the algorithms. Each parameter created will be based
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
        dbscan_param_pairs = dbscanTools.dbscan_param_space_search(self.parameters["evaluation"]["maximum_noise"], 
                                                                   self.distance_matrix)
        for (minpts, eps) in dbscan_param_pairs:
            run_parameter = ParametersGenerator.get_base_parameters()
            run_parameter["minpts"] = minpts
            run_parameter["eps"] = eps
            run_parameters.append(run_parameter)
        return run_parameters, []
    
        
'''
Created on 06/02/2013

@author: victor
'''

import pyproclust.algorithms.gromos.parametersGeneration as gromosParametersGeneration
import pyproclust.algorithms.kmedoids.parametersGeneration as kmedoidsParametersGeneration
import pyproclust.algorithms.random.parametersGeneration as randomParametersGeneration
import pyproclust.algorithms.spectral.parametersGeneration as spectralParametersGeneration
import pyproclust.algorithms.hierarchical.parametersGeneration as hierarchicalParametersGeneration
import pyproclust.algorithms.dbscan.parametersGeneration as dbscanParametersGeneration

class AlgorithmRunParametersGenerator(object):

    def __init__(self, parameters, matrix_handler):
        """
        Class creator.
        
        @param parameters: Script parameters.
        
        @param distance_matrix: The distance matrix we are using.
        """
        self.distance_matrix = matrix_handler.distance_matrix
        self.parameters = parameters
    
    def get_parameters_for_type(self, algorithm_type):
        """
        Tries to get some parameters to do a good exploration using one algorithm.
        
        @param algorithm_type: A valid algorithm type for which to get the parameters.
        
        @return: A tuple consisting on the parameters we are going to use and, if the clusterings were
        generated to get the clusterings, then those clusterings are also returned for the sake of performance.
        In this case the algorithm parameters in position 'i' correspond to the created cluster in position 'i'. 
        """
        generator = None
        
        if algorithm_type == "gromos":
            generator = gromosParametersGeneration.ParametersGenerator(self.parameters, self.distance_matrix)
                
        elif algorithm_type == "spectral":
            generator = spectralParametersGeneration.ParametersGenerator(self.parameters, self.distance_matrix)
                
        elif algorithm_type == "kmedoids":
            generator = kmedoidsParametersGeneration.ParametersGenerator(self.parameters, self.distance_matrix)
            
        elif algorithm_type == "hierarchical":
            generator = hierarchicalParametersGeneration.ParametersGenerator(self.parameters, self.distance_matrix)
                 
        elif algorithm_type == "dbscan":
            generator = dbscanParametersGeneration.ParametersGenerator(self.parameters, self.distance_matrix)
                 
        elif algorithm_type == "random":
            generator = randomParametersGeneration.ParametersGenerator(self.parameters, self.distance_matrix)
            
        return generator.get_parameters()

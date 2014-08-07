"""
Created on 06/02/2013

@author: victor
"""

import pyproct.clustering.algorithms.gromos.parametersGeneration as gromosParametersGeneration
import pyproct.clustering.algorithms.kmedoids.parametersGeneration as kmedoidsParametersGeneration
import pyproct.clustering.algorithms.random.parametersGeneration as randomParametersGeneration
import pyproct.clustering.algorithms.spectral.parametersGeneration as spectralParametersGeneration
import pyproct.clustering.algorithms.hierarchical.parametersGeneration as hierarchicalParametersGeneration
import pyproct.clustering.algorithms.dbscan.parametersGeneration as dbscanParametersGeneration

class AlgorithmRunParametersGenerator(object):

    def __init__(self, parameters, matrix_handler):
        """
        Class creator.
        
        @param parameters: Script parameters.
        
        @param distance_matrix: The distance matrix we are using.
        """
        self.matrix_handler = matrix_handler
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
            generator = gromosParametersGeneration.ParametersGenerator(self.parameters, self.matrix_handler)
                
        elif algorithm_type == "spectral":
            generator = spectralParametersGeneration.ParametersGenerator(self.parameters, self.matrix_handler)
                
        elif algorithm_type == "kmedoids":
            generator = kmedoidsParametersGeneration.ParametersGenerator(self.parameters, self.matrix_handler)
            
        elif algorithm_type == "hierarchical":
            generator = hierarchicalParametersGeneration.ParametersGenerator(self.parameters, self.matrix_handler)
                 
        elif algorithm_type == "dbscan":
            generator = dbscanParametersGeneration.ParametersGenerator(self.parameters, self.matrix_handler)
                 
        elif algorithm_type == "random":
            generator = randomParametersGeneration.ParametersGenerator(self.parameters, self.matrix_handler)
            
        return generator.get_parameters()

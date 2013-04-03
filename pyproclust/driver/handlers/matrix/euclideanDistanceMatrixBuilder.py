'''
Created on 13/02/2013

@author: victor
'''
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyRMSD.condensedMatrix import CondensedMatrix
import scipy.spatial.distance

class EuclideanDistanceMatrixBuilder(object):

    def __init__(self):
        pass
    
    @classmethod
    def build(cls, trajectory_handler, matrix_creation_parameters):
        """
        
        @param trajectory_handler:  
        @param matrix_creation_parameters: 
        
        @return: The created distances matrix.
        """
        #pdb = trajectory_handler.getJoinedPDB()
        
        # Build calculator with fitting coordinate sets
        fit_selection_coordsets = trajectory_handler.getSelection(matrix_creation_parameters["fit_selection"])
        calculator = RMSDCalculator(
                                    coordsets = fit_selection_coordsets,
                                    calculatorType = "QTRFIT_OMP_CALCULATOR", 
                                    modifyCoordinates = False)
        
        # Superpose iteratively
        calculator.iterativeSuperposition()
        
        #Then calculate distances
        body_selection_string = matrix_creation_parameters["body_selection"]
        body_selection_coordsets = trajectory_handler.getSelection(body_selection_string)
        
        trajectory_handler.setWorkingCoordinates(matrix_creation_parameters["body_selection_string"])
        
        return EuclideanDistanceMatrixBuilder(body_selection_coordsets)
    
    @classmethod
    def calculate_geom_center(cls, coordinates):
        """
        
        """
        # Calculate geom center
        centers = coordinates.mean(1)
        distances = scipy.spatial.distance.pdist(centers, 'euclidean')
        return  CondensedMatrix(distances)
        
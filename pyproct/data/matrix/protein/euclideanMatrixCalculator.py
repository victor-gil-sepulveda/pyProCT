"""
Created on 13/02/2013

@author: victor
"""
import scipy.spatial.distance
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyRMSD.condensedMatrix import CondensedMatrix

class EuclideanMatrixCalculator(object):
    
    CALCULATION_METHOD = "distance"
    
    def __init__(self):
        pass

    @classmethod
    def calculate(cls, trajectory_handler, matrix_parameters):
        """
        Will generate the CondensedMatrix filled with the all vs all geometric center distances of the "body_selection"
        coordinates (which will usually be a ligand).

        @param trajectory_handler: The handler containing selection strings, pdb info and coordsets.
        @param matrix_parameters: The creation parameters (from the initial script).
        
         - "distance": Euclidean distance of the geometrical center of one body.
    
                "parameters":{
                    "fit_selection":  String,
                    "body_selection": String,
                }
    
                "fit_selection": The Prody selection string used to describe the atoms to be superposed.
                "body_selection": Another Prody selection string that describes the element that will be used
                to get the euclidean distances.
        

        @return: The created distance matrix.
        """

        # Build calculator with fitting coordinate sets ...
        fit_selection_coordsets = trajectory_handler.getSelection(trajectory_handler.fitting_selection)

        # and calculation coordsets (we want them to be moved along with the fitting ones)
        body_selection_coordsets = trajectory_handler.getSelection(trajectory_handler.calculation_selection)

        calculator = RMSDCalculator(calculatorType = "QTRFIT_OMP_CALCULATOR",
                 fittingCoordsets = fit_selection_coordsets,
                 calculationCoordsets = body_selection_coordsets)

        # Superpose iteratively (will modify all coordinates)
        calculator.iterativeSuperposition()

        # Working coordinates are changed to the body coordinates (to be used later for instance
        # with clustering metrics)
        trajectory_handler.setWorkingCoordinates(trajectory_handler.calculation_selection)
        distances = cls.calculate_geom_center(body_selection_coordsets)
        matrix = CondensedMatrix(distances)
        return matrix

    @classmethod
    def calculate_geom_center(cls, coordinates):
        """
        Generates a condensed matrix with the euclidean distances between the geometrical centers of the conformations passed as
        input.

        @param coordinates: Coordinates set from which calculating the geometrical centers (one geometrical center per conformation).

        @return: The contents of the condensed matrix resulting of calculating all euclidean distances between the aforemetioned centers.
        """
        # Calculate geom centers
        centers = coordinates.mean(1)
        distances = scipy.spatial.distance.pdist(centers, 'euclidean')
        return  distances

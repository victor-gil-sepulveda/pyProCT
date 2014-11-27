'''
Created on 27/11/2014

@author: victor
'''
from pyRMSD.RMSDCalculator import RMSDCalculator

class euclideanDistanceBuilder(object):

    def __init__(self, params):
        pass
    
    @classmethod
    def build(cls, data_handler, matrix_params):
        """
        - "distance": Euclidean distance of the geometrical center of one body.
    
                "parameters":{
                    "fit_selection":  String,
                    "body_selection": String,
                }
    
                "fit_selection": The Prody selection string used to describe the atoms to be superposed.
                "body_selection": Another Prody selection string that describes the element that will be used
                to get the euclidean distances.
        
        """
        # Build calculator with fitting coordinate sets ...
        fit_selection_coordsets = data_handler.get_data().getFittingCoordinates()

        # and calculation coordsets (we want them to be moved along with the fitting ones)
        body_selection_coordsets = data_handler.get_data().getCalculationCoordinates()

        calculator = RMSDCalculator(calculatorType = "QTRFIT_OMP_CALCULATOR",
                 fittingCoordsets = fit_selection_coordsets,
                 calculationCoordsets = body_selection_coordsets)

        # Superpose iteratively (will modify all coordinates)
        calculator.iterativeSuperposition()
        
        centers = body_selection_coordsets.mean(1)
        
        return centers

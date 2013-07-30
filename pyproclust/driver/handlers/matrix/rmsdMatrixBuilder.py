'''
Created on 13/02/2013

@author: victor
'''
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyRMSD.condensedMatrix import CondensedMatrix

class RMSDMatrixBuilder(object):

    def __init__(self):
        pass
    
    @classmethod
    def build(cls, trajectory_handler, matrix_creation_parameters):
        """
        Generates a matrix with the method used in the handler creation.
        
        @param trajectory_handler: 
        @param matrix_creation_parameters: 
        
        @return: The created matrix.
        """

        fit_selection_string = matrix_creation_parameters["fit_selection"]
        fit_selection_coordsets = trajectory_handler.getSelection(fit_selection_string)
        trajectory_handler.setWorkingCoordinates(fit_selection_string)


        calculator = RMSDCalculator(calculatorType = "QTRFIT_OMP_CALCULATOR",
                                    fittingCoordsets  = fit_selection_coordsets)
        
        # Apply calculation selection if needed
        calc_selection_string = matrix_creation_parameters["calc_selection"]
        if calc_selection_string != "" and calc_selection_string != fit_selection_string:
            calc_selection_coordsets = trajectory_handler.getSelection(calc_selection_string)
            trajectory_handler.setWorkingCoordinates(calc_selection_string)
            calculator.setCalculationCoordinates(calc_selection_coordsets)
            calculator = RMSDCalculator(calculatorType = "QTRFIT_OMP_CALCULATOR",
                                        fittingCoordsets  = fit_selection_coordsets,
                                        calculationCoordsets = calc_selection_coordsets)
        
        rmsds = calculator.pairwiseRMSDMatrix()
        
        return CondensedMatrix(rmsds)
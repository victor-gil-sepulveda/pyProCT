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

        fit_selection_coordsets = trajectory_handler.getSelection(matrix_creation_parameters["fit_selection"])
        trajectory_handler.setWorkingCoordinates(matrix_creation_parameters["fit_selection"])
        calculator = RMSDCalculator(coordsets = fit_selection_coordsets,
                                    calculatorType = "QTRFIT_OMP_CALCULATOR", 
                                    modifyCoordinates = False)
        
        # Apply calculation selection if needed
        calc_selection_string = matrix_creation_parameters["calc_selection"]
        if calc_selection_string != "":
            calc_selection_coordsets = trajectory_handler.getSelection(matrix_creation_parameters["calc_selection"])
            trajectory_handler.setWorkingCoordinates(calc_selection_string)
            calculator.setCalculationCoordinates(calc_selection_coordsets)
        
        rmsds = calculator.pairwiseRMSDMatrix()
        
        return CondensedMatrix(rmsds)
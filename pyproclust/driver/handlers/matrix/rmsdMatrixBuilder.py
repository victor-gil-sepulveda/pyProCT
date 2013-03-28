'''
Created on 13/02/2013

@author: victor
'''
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyRMSD.condensedMatrix import ConsensedMatrix

class RMSDMatrixBuilder(object):

    def __init__(self):
        pass
    
    def build(self, trajectory_handler, matrix_creation_parameters):
        """
        Generates a matrix with the method used in the handler creation.
        
        @param trajectory_handler: 
        @param matrix_creation_parameters: 
        
        @return: The created matrix.
        """
        pdb = trajectory_handler.getJoinedPDB()
        
        # Build calculator with fitting coordinate sets
        fit_selection_string = matrix_creation_parameters["fit_selection"]
        fit_selection_coordsets = None
        if fit_selection_string == "":
            fit_selection_coordsets = pdb.getCoordsets()
        else:
            fit_selection_coordsets = pdb.select(fit_selection_string).
        calculator = RMSDCalculator(    
                                    coordsets = fit_selection_coordsets, 
                                    calculatorType = "QTRFIT_OMP_CALCULATOR", 
                                    modifyCoordinates = False)
        
        # Apply calculation selection if needed
        calc_selection_string = matrix_creation_parameters["calc_selection"]
        if calc_selection_string != "":
            calc_selection_coordsets = pdb.select(calc_selection_string)
            calculator.setCalculationCoordinates(calc_selection_coordsets)
        
        rmsds = calculator.pairwiseRMSDMatrix()
        
        return CondensedMatrix(rmsds)
"""
Created on Sep 3, 2014

@author: victor
"""

class ProteinEnsembleData(object):
    """
    Holds a protein ensemble and helps to handle it (selections etc...)
    """
    
    def __init__(self,  structure_ensemble, selection_params):
        self.structure_ensemble = structure_ensemble
        self.handle_selection_parameters(selection_params)
    
    def getSelection(self, selection_string):
        """
        Returns a Prody modifiable selection.
        """
        return self.structure_ensemble.select(selection_string).copy()
    
    def getSelectionCoordinates(self, selection_string):
        """
        Returns the coordinates of an arbitrary selection.
        """
        return self.getSelection(selection_string).getCoordsets()

    def getFittingCoordinates(self):
        """
        Returns the coordinates for the fitting selection.
        """
        if self.fitting_selection is not None:
            return self.getSelectionCoordinates(self.fitting_selection)
        else:
            return None

    def getCalculationCoordinates(self):
        """
        Returns the coordinates for the calculation selection.
        """
        if self.calculation_selection is not None:
            return self.getSelectionCoordinates(self.calculation_selection)
        else:
            return None
    
    def handle_selection_parameters(self, selection_parameters):
        """
        Helper funtion to handle selection parameters (different parameter names can have almost the same
        functionality and are treated internally in the same way).

        @param matrix_parameters: The parameters chunk that controls matrix selections.
        """
        # Store the main selections we can do
        self.fitting_selection = self.calculation_selection = None

        if "fit_selection" in selection_parameters:
            self.fitting_selection = selection_parameters["fit_selection"]

        if "dist_fit_selection" in selection_parameters:
            self.fitting_selection = selection_parameters["dist_fit_selection"]

        if "calc_selection" in selection_parameters:
            self.calculation_selection = selection_parameters["calc_selection"]

        if "body_selection" in selection_parameters:
            self.calculation_selection = selection_parameters["body_selection"]
        
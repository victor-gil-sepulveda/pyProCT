"""
Created on Sep 3, 2014

@author: victor
"""
from pyproct.data.handler.data import Data
from pyproct.tools.prodyTools import removeAllCoordsetsFromStructure
from prody.measure.measure import calcPhi, calcPsi

class ProteinEnsembleData(Data):
    """
    Holds a protein ensemble with ALL conformations and helps to handle it (selections etc...).
    TODO: Methods to lowercase.
    """
    
    def __init__(self,  structure_ensemble, model_numbers, model_remarks, selection_params):
        # Add number of elements to data instead of handler
        self.structure_ensemble = structure_ensemble
        self.handle_selection_parameters(selection_params)
        self.model_numbers = model_numbers
        self.model_remarks = model_remarks

    def get_element(self, element_id):
        """
        We must override this guy. In this case the behaviour is to obtain a 
        full Prody structure.
        """
        element_coordinates = self.structure_ensemble.getCoordsets()[element_id]
        structure = self.structure_ensemble.copy()
        removeAllCoordsetsFromStructure(structure)
        structure.addCoordset(element_coordinates)
        return structure
    
    def get_elements(self, elements_list):
        """
        When we work with ensembles it is better to return a single prody ensemble with all
        the elements we want.
        Returns a copy! Modifying it does not modify the stored structure.
        """
        element_coordinates = self.structure_ensemble.getCoordsets()[elements_list]
        structure = self.structure_ensemble.copy()
        removeAllCoordsetsFromStructure(structure)
        for coordset in element_coordinates:
            structure.addCoordset(coordset)
        return structure
    
    def get_all_elements(self):
        """
        """
        return self.structure_ensemble
    
    def get_number_of_elements(self):
        return self.structure_ensemble.numCoordsets()
    
    def getSelection(self, selection_string):
        """
        Returns a Prody modifiable selection (is a copy).
        """
        return self.structure_ensemble.select(selection_string).copy()
    
    def getCoordinates(self):
        """
        Returns all the coordinates of the ensemble.
        """
        return self.getSelection("all").getCoordsets()
    
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

        @param selection_parameters: The parameters chunk that controls matrix selections.
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
        
        # If both selections are the sa
        if self.fitting_selection == self.calculation_selection:
            self.calculation_selection = None 
        
    def get_dihedrals_for_conformation(self, conformation):
        coordsets = self.structure_ensemble.getCoordsets()
        
        self.structure_ensemble.setCoords(coordsets[conformation])
        
        dihedral_angles = []
        for residue in self.structure_ensemble.iterResidues():
            try:
                dihedral_angles.append(calcPhi(residue, radian=False, dist=None))
            except ValueError:
                dihedral_angles.append(0)

            try:
                dihedral_angles.append(calcPsi(residue, radian=False, dist=None))
            except ValueError:
                dihedral_angles.append(0)

        # 0 links with Nth residue and Nth with 0th. Those values are not needed anyway.
        return dihedral_angles[1:-1]
    
    def get_all_remarks(self):
        return self.model_remarks
    
    def get_all_model_numbers(self):
        return self.model_numbers

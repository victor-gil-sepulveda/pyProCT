"""
Created on 19/09/2012

@author: victor
"""
import pyproct.tools.commonTools as common
from pyproct.driver.observer.observable import Observable
import prody
import os.path
from pyproct.tools.prodyTools import removeAllCoordsetsFromStructure
from pyproct.tools.commonTools import get_parameter_value

class ProteinStructureEnsembleDataLoader(Observable):
    LOADER_TYPE = "protein::ensemble"
    
    def __init__(self, data_parameters, observer):
        """
        """
        self.structures = []
        self.structure_ensemble = None
        self.handle_selection_parameters(get_parameter_value("matrix", {"fit_selection":"all"}))

    def load(self, data_source):
        """
        Adds a new ensemble.
        """
        self.structures.append(self.load_structure_ensemble(data_source))
    
    def close(self):
        """
        Prepares the merged structure.
        """
        if len(self.structures) == 0:
            common.print_and_flush("[ERROR ProteinStructureEnsembleData:close] No loaded structures. Exiting...\n")
            exit()

        self.structure_ensemble = self.generate_merged_structure_ensemble(self.structures) 
        
        # We free some memory
        del self.structures
        
        print "%d conformations of %d atoms were read."% (self.structure_ensemble.numCoordsets(),self.structure_ensemble.numAtoms())

    def load_structure_ensemble(self, source):
        """
        Loads a structure file (pdb or dcd) and updates source info.

        :param source: Is a DataSource object with one of this sets of keywords:
        
        - For 'pdb' files:

            {
                "source": ... ,
                "base_selection": ...
            }

        Where 'file' contains the path of the pdb file we want to load.

        - For 'dcd' files:

            {
                "source": ...,
                "atoms_source": ...,
                "base_selection": ...
            }

        Where 'file' contains the path of the 'dcd' file we want to load and 'atoms_file' the source of the pdb file containing
        the atomic information.

        In both cases 'base_selection' is a Prody selection string that performs an initial selection of the atoms. This is
        useful when we want to load more than one file with different number of atoms and its goal is to allow the selection
        of the common atoms. It is up to the user to maintain a 1 to 1 mapping between the atoms of each of the files.

        The source object will be enriched with some extra information from the loaded structure ensemble.

        :return: Prody's structure object with the loaded ensemble
        """
        _, ext = os.path.splitext(source.get_path())
        
        if ext == ".dcd":
            structure = prody.parsePDB(source.get_info("atoms_source"))
            # Leave only atomic information
            removeAllCoordsetsFromStructure(structure)
            dcd_data = prody.DCDFile(source.get_path())
            coordsets = dcd_data.getCoordsets()
            # Add all coordsets to atomic information
            for coordset in coordsets:
                structure.addCoordset(coordset)
            
        elif ext == ".pdb":
            structure = prody.parsePDB(path)
        else:
            print "[ERROR][ProteinStructureEnsembleData::get_structure] pyProCT does not know hot to load the file %s (unknown extension '%s')"%(path,ext)
            exit()
        
        if source.has_info("base_selection"):
            structure = structure.select(source.get_info("base_selection")).copy()
            if structure is None:
                common.print_and_flush("[ERROR ProteinStructureEnsembleData::get_structure] Improductive base selection (%s). Exiting...\n"%source.get_info("base_selection"))
                exit()

        source.add_info("number of conformations", structure.numCoordsets())
        source.add_info(["number of atoms"], structure.numAtoms())
        
        return  structure

    def generate_merged_structure_ensemble(self, structures):
        """
        Merges all handled structures into a single Prody AtomGroup object.

        @return: The prody object with all read coordsets for certain selection.
        """
        merged_ensemble = None
        for structure in structures:
            if merged_ensemble is None:
                merged_ensemble = structure
            else:
                for coordset in structure.getCoordsets():
                    merged_ensemble.addCoordset(coordset)
        return merged_ensemble

    def getSelection(self, selection_string):
        """
        """
        return self.structure_ensemble.select(selection_string).copy().getCoordsets()

    def getFittingCoordinates(self):
        if self.fitting_selection is not None:
            return self.getSelection(self.fitting_selection)
        else:
            return None

    def getCalculationCoordinates(self):
        if self.calculation_selection is not None:
            return self.structure_ensemble(self.calculation_selection)
        else:
            return None
    
    def handle_selection_parameters(self, matrix_parameters):
        """
        Helper funtion to handle selection parameters (different parameter names can have almost the same
        functionality and are treated internally in the same way).

        @param matrix_parameters: The parameters chunk that controls matrix selections.
        """
        # Store the main selections we can do
        self.fitting_selection = self.calculation_selection = None

        if "fit_selection" in matrix_parameters:
            self.fitting_selection = matrix_parameters["fit_selection"]

        if "dist_fit_selection" in matrix_parameters:
            self.fitting_selection = matrix_parameters["dist_fit_selection"]

        if "calc_selection" in matrix_parameters:
            self.calculation_selection = matrix_parameters["calc_selection"]

        if "body_selection" in matrix_parameters:
            self.calculation_selection = matrix_parameters["body_selection"]
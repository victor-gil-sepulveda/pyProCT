"""
Created on 19/09/2012

@author: victor
"""
import pyproct.tools.commonTools as common
import pyproct.tools.pdbTools as pdb_tools
import prody
import os.path
from pyproct.tools.prodyTools import removeAllCoordsetsFromStructure
from pyproct.data.handler.protein.proteinEnsembleData import ProteinEnsembleData
from pyproct.data.handler.dataLoader import DataLoader

class ProteinEnsembleDataLoader(DataLoader):
    """
    A loader for pdb and dcd ensemble files. As a DataLoader, it must implement
    its own load and close methods 
    """
    
    LOADER_TYPE = "protein::ensemble"
    
    def __init__(self, data_params):
        super(ProteinEnsembleDataLoader, self).__init__(data_params)
        self.model_numbers = []
        self.model_remarks = []
        
    def close(self):
        """
        Prepares the merged structure and returns the ensemble data object
        """
        super(ProteinEnsembleDataLoader, self).close()
        
        structure_ensemble = self.generate_merged_structure_ensemble(self.loaded_data) 
        
        # We free some memory
        del self.loaded_data
        
        print "%d conformations of %d atoms were read."% (structure_ensemble.numCoordsets(),
                                                          structure_ensemble.numAtoms())
        
        return ProteinEnsembleData(structure_ensemble,
                                   self.model_numbers,
                                   self.model_remarks, 
                                   self.data_params.get_value("matrix.parameters", 
                                                       {"fit_selection": "all"}))

    def generate_merged_structure_ensemble(self, structures):
        """
        Merges all handled structures into a single Prody AtomGroup object.

        @return: The prody object with all read coordsets for certain selection.
        
        TODO: If atoms of the different loaded files are not ordered, then selections
        (which refer to the first loaded model) will not work as expected.
        """
        merged_ensemble = None
        for structure in structures:
            if merged_ensemble is None:
                merged_ensemble = structure
            else:
                for coordset in structure.getCoordsets():
                    merged_ensemble.addCoordset(coordset)
        return merged_ensemble

    def load_data_from_source(self, source):
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
            structure = prody.parsePDB(source.get_path())
        else:
            print "[ERROR][ProteinStructureEnsembleData::get_structure] pyProCT does not know hot to load the file %s (unknown extension '%s')"%(source.get_path(),ext)
            exit()
        
        if source.has_info("base_selection"):
            structure = structure.select(source.get_info("base_selection")).copy()
            if structure is None:
                common.print_and_flush("[ERROR ProteinStructureEnsembleData::get_structure] Improductive base selection (%s). Exiting...\n"%source.get_info("base_selection"))
                exit()

        source.add_info("number_of_conformations", structure.numCoordsets())
        source.add_info("number_of_atoms", structure.numAtoms())
        
        self.model_numbers.extend(self.get_model_numbers(source, structure.numCoordsets()))
        self.model_remarks.extend(self.get_remarks(source, structure.numCoordsets()))
        
        return  structure, structure.numCoordsets()


    def get_model_numbers(self, source, number_of_conformations):
        """
        :param source:  
        """
        _, ext = os.path.splitext(source.get_path())
        
        if ext == ".dcd":
            return range(1, number_of_conformations+1)
        else:
            model_lines = pdb_tools.get_model_tags(source.get_path())
            if len(model_lines) != number_of_conformations:
                # by default
                return range(1, number_of_conformations+1)
            else:
                return [int(model.split()[1]) for model in model_lines]
            
    def get_remarks(self, source, number_of_conformations):
        _, ext = os.path.splitext(source.get_path())
        
        if ext == ".dcd":
            return [[]]*number_of_conformations
        else:
            remark_groups =  pdb_tools.get_remarks(source.get_path())
            if len(remark_groups) != number_of_conformations:
                # TODO: WARNING HERE
                return [[]]*number_of_conformations
            else:
                return remark_groups
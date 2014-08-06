"""
Created on 19/09/2012

@author: victor
"""
import pyproct.tools.commonTools as common
from pyproct.driver.observer.observable import Observable
import prody
import os.path
from pyproct.tools.prodyTools import removeAllCoordsetsFromStructureLeavingFirst
from pyproct.driver.parameters import ProtocolParameters

class TrajectoryHandler(Observable):

    def __init__(self, parameters, observer):
        """
        Class creator. It parses the needed files and extracts info and coordinates.
        """

        super(TrajectoryHandler,self).__init__(observer)

        print "Reading conformations..."
        prody.confProDy(verbosity="none")

        self.parameters = parameters
        matrix_parameters = parameters.get_value("data.matrix.parameters", default_value=ProtocolParameters.empty())
        parameters["data"]["files"] = self.expand_file_lists(parameters["data"]["files"])
        self.files = parameters["data"]["files"]
        self.pdbs = []

        if len(self.files) == 0:
            common.print_and_flush( "[ERROR] no pdbs. Exiting...\n")
            self.notify("SHUTDOWN","No pdbs defined in script.")
            exit()

        self.notify("Loading","Loading Trajectories")

        # Bookmarking structure
        self.bookmarking = {
                             "pdb": None,
                             "selections": {}
        }

        merged_structure = self.getMergedStructure()
        self.coordsets = merged_structure.getCoordsets()
        self.number_of_conformations = self.coordsets.shape[0]
        self.number_of_atoms = self.coordsets.shape[1]

        self.handle_selection_parameters(matrix_parameters)
        print "%d conformations of %d atoms were read."% (merged_structure.numCoordsets(),merged_structure.numAtoms())

    def expand_file_lists(self, files):
        """
        Loads the PDB files specified into a file (maintaining order).
        """
        tmp_file_list = []
        for file_info in files:
            if isinstance(file_info, basestring):
                # Then it must be a txt file
                name, ext = os.path.splitext(file_info)
                if ext in [".txt", ".list"]:
                    # Load the file
                    lines = open(file_info,"r").readlines()
                    for line in lines:
                        l = line.strip()
                        if "," in l:
                            file_name, selection = l.split(",")
                            tmp_file_list.append({"file":file_name,"base_selection":selection})
                        else:
                            tmp_file_list.append(l)
                else:
                    tmp_file_list.append(file_info)
            else:
                tmp_file_list.append(file_info)
        return tmp_file_list


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

    def check_extension(self, ext):
        """
        Helper function to check if the file extension is allowed. If not it shuts down the program.

        @param ext: The extension string (with the separating period!)

        @return: Nothing (exits if the condition is not fulfilled)
        """
        if not ext in [".dcd",".pdb"]:
            common.print_and_flush( "[ERROR] pyProCT cannot read this file format.\n")
            self.notify("SHUTDOWN","Wrong file format.")
            exit()

    def get_structure(self, file_info):
        """
        Loads a structure file (pdb or dcd) and fills its structure_info data for logging.

        @param file_info: Is a string containing the path of the file or a dictionary with this structure:
        'pdb' files:

            {
                "file": ... ,
                "base_selection": ...
            }

        Where 'file' contains the path of the pdb file we want to load.

        'dcd' files:

            {
                "file": ...,
                "atoms_file": ...,
                "base_selection": ...
            }

        Where 'file' contains the path of the 'dcd' file we want to load and 'atoms_file' the source of the pdb file containing
        the atomic information.

        In both cases 'base_selection' is a Prody selection string that performs an initial selection of the atoms. This is
        useful when we want to load more than one file with different number of atoms and its goal is to allow the selection
        of the common atoms. It is up to the user to maintain a 1 to 1 mapping between the atoms of each of the files.

        @return: A tuple containing the structure object and a structure_info dictionary.
        """
        structure_info = {
              "source":"",
              "source of atoms":"",
              "base selection": "",
              "number of conformations": "",
              "number of atoms":  ""
        }

        if isinstance(file_info, basestring):
            # Then is a path, and must be a pdb
            path = file_info
            structure_info["source"] = path

            name, ext = os.path.splitext(path)

            self.check_extension(ext)

            if ext == ".dcd":
                common.print_and_flush( "[ERROR TrajectoryHandler::get_structure] Path format can only be used with pdb files. Exiting...\n")
                self.notify("SHUTDOWN", "Fatal error reading pdbs.")
                exit()
            else:
                structure = prody.parsePDB(path)
                structure_info["number of conformations"] = structure.numCoordsets()
                structure_info["number of atoms"] = structure.numAtoms()
                return  structure, structure_info
        else:
            # {"file":  , "selection":  } object or
            # {"file": , "atoms_file":, "selection"} if the file is a dcd file
            path = file_info["file"]
            structure_info["source"] = path
            name, ext = os.path.splitext(path)
            self.check_extension(ext)

            if ext == ".dcd":
                structure_info["source of atoms"] = file_info["atoms_file"]

                structure = prody.parsePDB(file_info["atoms_file"])
                removeAllCoordsetsFromStructureLeavingFirst(structure)
                dcd_data = prody.DCDFile(path)
                coordsets = dcd_data.getCoordsets()

                for coordset in coordsets:
                    structure.addCoordset(coordset)
            else:
                structure = prody.parsePDB(path)

            if "base_selection" in file_info and file_info["base_selection"] !=  "":
                structure = structure.select(file_info["base_selection"])
                structure_info["base selection"]=file_info["base_selection"]

            structure_info["number of conformations"] = structure.numCoordsets()
            structure_info["number of atoms"] = structure.numAtoms()
            return  structure, structure_info

    def getMergedStructure(self):
        """
        Merges all handled structures into a single Prody AtomGroup object. If there's any error, the program must exit, and
        any thread must be stopped.

        @return: The prody object with all read coordsets for certain selection.
        """
        merged_pdb = None
        if self.bookmarking["pdb"] is None:
            try:
                for file_info in self.files:
                    structure, structure_info = self.get_structure(file_info)
                    self.pdbs.append(structure_info)

                    if merged_pdb is None:
                        merged_pdb = structure.copy()
                    else:
                        for coordset in structure.getCoordsets():
                            merged_pdb.addCoordset(coordset)
                self.bookmarking["pdb"] = merged_pdb
            except Exception, e:
                print "[ERROR TrajectroyHandler::getMergedStructure] fatal error reading pdbs.\nError: %s\n Program will halt now ..."%e.message
                self.notify("SHUTDOWN", "Fatal error reading pdbs.")
                exit()

        return self.bookmarking["pdb"]

    def setWorkingCoordinates(self, selection_string):
        self.bookmarking["working"] = selection_string

    def getWorkingCoordinates(self):
        selection_string = self.bookmarking["working"]

        if selection_string == "":
            return self.getMergedStructure().getCoordsets()

        if not selection_string in self.bookmarking["selections"]:
            return self.getSelection(selection_string)
        else:
            return self.bookmarking["selections"][selection_string];

    def getSelection(self, selection_string):
        if self.bookmarking["pdb"] is None:
            self.getMergedStructure()

        pdb = self.bookmarking["pdb"]

        if not selection_string in self.bookmarking["selections"]:
            selection_coordsets = None
            if selection_string == "":
                selection_coordsets = pdb.getCoordsets()
            else:
                selection_coordsets = pdb.select(selection_string).getCoordsets()

            self.bookmarking["selections"][selection_string] = selection_coordsets

        else:
            selection_coordsets = self.bookmarking["selections"][selection_string]

        return selection_coordsets

    def getFittingCoordinates(self):
        if self.fitting_selection is not None:
            return self.getSelection(self.fitting_selection)
        else:
            return None

    def getCalculationCoordinates(self):
        if self.calculation_selection is not None:
            return self.getSelection(self.calculation_selection)
        else:
            return None
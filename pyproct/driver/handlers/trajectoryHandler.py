'''
Created on 19/09/2012

@author: victor
'''
import pyproct.tools.commonTools as common
from pyproct.driver.observer.observable import Observable
import prody
import os.path

class TrajectoryHandler(Observable):

    def __init__(self, parameters, observer):
        """
        Class creator. It parses the needed files and extracts info and coordinates.
        """
        super(TrajectoryHandler,self).__init__(observer)

        matrix_parameters = parameters["data"]["matrix"]['parameters']
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

        self.coordsets = self.getMergedPDB().getCoordsets()
        self.number_of_conformations = self.coordsets.shape[0]
        self.number_of_atoms = self.coordsets.shape[1]

        self.handle_selection_paramters(matrix_parameters)


    def handle_selection_paramters(self, matrix_parameters):
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

    @classmethod
    def extract_file_data(cls, path, structure):
        """
        Creates a pdb dictionary with source, number of frames and number of atoms.

        @param pdb: The full path of the pdb from which to extract data.

        @return: The aforementioned dictionary.
        """
        return {
                  "source":path,
                  "conformations": structure.numCoordsets(),
                  "atoms":  structure.numAtoms()
        }

    def get_structure_file(self, path):
        name, ext = os.path.splitext(path)

        if not ext in [".dcd",".pdb"]:
            common.print_and_flush( "[ERROR] pyProCT cannot read this file format.\n")
            self.notify("SHUTDOWN","Wrong file format.")
            exit()

        if ext == ".dcd":
            return prody.DCDFile(path)
        else:
            return  prody.parsePDB(path)

    def getMergedPDB(self):
        """
        Merges all handled pdbs into a single Prody pdb object. If there's any error, the program must exit, and
        any thread must be stopped.

        @return: The prody object with all read coordsets for certain selection.
        """
        merged_pdb = None
        if self.bookmarking["pdb"] is None:
            try:
                for path in self.files:

                    pdb = self.get_structure_file(path)
                    self.pdbs.append(self.extract_file_data(path, pdb))

                    if merged_pdb is None:
                        merged_pdb = pdb
                    else:
                        for coordset in pdb.getCoordsets():
                            merged_pdb.addCoordset(coordset)
                self.bookmarking["pdb"] = merged_pdb
            except Exception, e:
                print "[ERROR TrajectroyHandler::getMergedPDB] fatal error reading pdbs.\nError: %s\n Program will halt now ..."%e.message
                self.notify("SHUTDOWN", "Fatal error reading pdbs.")
                exit()

        return self.bookmarking["pdb"]

    def setWorkingCoordinates(self, selection_string):
        self.bookmarking["working"] = selection_string

    def getWorkingCoordinates(self):
        selection_string = self.bookmarking["working"]

        if selection_string == "":
            return self.getMergedPDB().getCoordsets()

        if not selection_string in self.bookmarking["selections"]:
            return self.getSelection(selection_string)
        else:
            return self.bookmarking["selections"][selection_string];

    def getSelection(self, selection_string):
        if self.bookmarking["pdb"] is None:
            self.getMergedPDB()

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
'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
import pyproclust.tools.pdbTools as pdb_tools
from pyproclust.driver.observer.observable import Observable
import prody

class TrajectoryHandler(Observable):

    def __init__(self, global_parameters, matrix_parameters, observer):
        """
        """
        super(TrajectoryHandler,self).__init__(observer)
        
        if len(global_parameters["pdbs"]) > 0:
            self.pdbs = [self.get_pdb_data(pdb_path) for pdb_path in global_parameters["pdbs"]]
        else:
            common.print_and_flush( "[ERROR] no pdbs. Exiting...\n")
            self.notify("SHUTDOWN","No pdbs defined in script.")
            exit()
            
        self.notify("Loading","Loading Trajectories")
        
        # Bookmarking
        self.bookmarking = { 
                             "pdb": None,
                             "selections": {}
        }
    
        self.coordsets = self.getJoinedPDB().getCoordsets()
        self.number_of_conformations = self.coordsets.shape[0]
        self.number_of_atoms = self.coordsets.shape[1]
        
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
    def get_pdb_data(cls, pdb):
        """
        Creates a pdb dictionary with source, number of frames and number of atoms.
        
        @param pdb: The full path of the pdb from which to extract data.
        
        @return: The aforementioned dictionary.
        """
        return {
                  "source":pdb,
                  "conformations": pdb_tools.get_number_of_frames(pdb),
                  "atoms":  pdb_tools.get_number_of_atoms(pdb)
        }
    
    def getJoinedPDB(self):
        """
        Merges all handled pdbs into a single Prody pdb object. If there's any error, the program must exit, and
        any thread must be stopped.
        
        @return: The prody object with all read coordsets for certain selection.
        """
        merged_pdb = None
        if self.bookmarking["pdb"] is None:
            try:
                for pdb_data in self.pdbs:
                    pdb = prody.parsePDB(pdb_data["source"])
                    if merged_pdb is None:
                        merged_pdb = pdb
                    else:
                        for coordset in pdb.getCoordsets():
                            merged_pdb.addCoordset(coordset)
                self.bookmarking["pdb"] = merged_pdb
            except Exception, e:
                print "[ERROR TrajectroyHandler::getJoinedPDB] fatal error reading pdbs.\nError: %s\n Program will halt now ..."%e.message
                self.notify("SHUTDOWN", "Fatal error reading pdbs.")
                exit()
        
        return self.bookmarking["pdb"]
    
    def setWorkingCoordinates(self, selection_string):
        self.bookmarking["working"] = selection_string
    
    def getWorkingCoordinates(self):
        selection_string = self.bookmarking["working"]
        
        if selection_string == "":
            return self.getJoinedPDB().getCoordsets()
        
        if not selection_string in self.bookmarking["selections"]:
            return self.getSelection(selection_string)
        else:
            return self.bookmarking["selections"][selection_string];
    
    def getSelection(self, selection_string):
        if self.bookmarking["pdb"] is None:
            self.getJoinedPDB()
        
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
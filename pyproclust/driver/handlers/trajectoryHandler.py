'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
import pyproclust.tools.pdbTools as pdb_tools
from pyRMSD.utils.proteinReading import Reader
from pyproclust.driver.observer.observable import Observable
import prody 

class TrajectoryHandler(Observable):

    def __init__(self, parameters, observer):
        super(TrajectoryHandler,self).__init__(observer)
        
        if parameters["global"]["action"] in ["comparison", "clustering"]:
            self.pdbs = [self.__do_pdb_dic(pdb_path) for pdb_path in parameters["global"]["pdbs"]]
        else:
            common.print_and_flush( "[ERROR] not known action. Exiting...\n")
            exit()
            
        reader = Reader("PRODY_READER")
        for pdb_description in self.pdbs:
            reader = reader.readThisFile(pdb_description["source"])
        
        self.notify("Loading","Loading Trajectories")
        
        self.coordsets = reader.read()
        
        self.number_of_conformations = reader.numberOfFrames
        self.number_of_atoms = reader.numberOfAtoms
        
    def __do_pdb_dic(self,pdb):
        return {
                  "source":pdb,
                  "conformations": pdb_tools.get_number_of_frames(pdb)
        }
        
    def getJoinedPDB(self):
        """
        Merges all handled pdbs into a single prody pdb object. If there's any error, the program must exit.
        
        @return: The prody object with all read coordsets.
        """
        merged_pdb = None
        try:
            for pdb_data in self.pdbs:
                pdb = prody.parsePDB(pdb_data["source"])
                if merged_pdb is None:
                    merged_pdb = pdb
                else:
                    merged_pdb.addCoordsets(pdb.getCoordsets())
        except Exception:
            print "[ERROR TrajectroyHandler::getJoinedPDB] fatal error reading pdbs. Program will exit..."
            self.notify("SHUTDOWN", "Fatal error reading pdbs.")
            exit()
        return merged_pdb
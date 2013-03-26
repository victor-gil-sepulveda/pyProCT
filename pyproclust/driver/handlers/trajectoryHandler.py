'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
import pyproclust.tools.pdbTools as pdb_tools
from pyRMSD.utils.proteinReading import Reader
from pyproclust.driver.observer.observable import Observable

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
'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
from pyRMSD.utils.proteinReading import Reader
from pyproclust.tools.pdbTools import get_number_of_frames
from pyproclust.protocol.observer.observable import Observable

class TrajectoryHandler(Observable):

    def __init__(self, parameters, observer):
        super(TrajectoryHandler,self).__init__(observer)
        
        if parameters["matrix"]["action"] == "clustering":
            self.pdbs = [self.__do_pdb_dic(parameters["matrix"]["pdb1"])]
            
        elif parameters["matrix"]["action"] == "comparison":
            self.pdbs = [self.__do_pdb_dic(parameters["matrix"]["pdb1"]),
                         self.__do_pdb_dic(parameters["matrix"]["pdb2"])]
            
        else:
            common.print_and_flush( "[ERROR] not known action. Exiting...\n")
            exit()
            
        reader = Reader("PRODY_READER")
        for pdb_description in self.pdbs:
            reader = reader.readThisFile(pdb_description["source"])
        
        self.notify("Loading","Loading Trajectories")
        if parameters["matrix"]["only_ca"]:
            reader = reader.gettingOnlyCAs()
        
        self.coordsets = reader.read()
        
        self.number_of_conformations = reader.numberOfFrames
        self.number_of_atoms = reader.numberOfAtoms
        
    def __do_pdb_dic(self,pdb):
        return {
                  "source":pdb,
                  "conformations": get_number_of_frames(pdb)
        }
'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
from pyproclust.tools.pdbTools import get_number_of_frames
from pyRMSD.utils.proteinReading import Reader

class TrajectoryHandler(object):

    def __init__(self, parameters):
        if parameters["matrix"]["action"] == "clustering":
            self.pdbs = [self.do_pdb_dic(parameters["matrix"]["pdb1"])]
            
        elif parameters["matrix"]["action"] == "comparison":
            self.pdbs = [self.do_pdb_dic(parameters["matrix"]["pdb1"]),
                         self.do_pdb_dic(parameters["matrix"]["pdb2"])]
            
        else:
            common.print_and_flush( "[ERROR] not known action. Exiting...\n")
            exit()
            
        reader = Reader("PRODY_READER")
        for pdb_description in self.pdbs:
            reader = reader.readThisFile(pdb_description["source"])
        
        if parameters["matrix"]["only_ca"]:
            reader = reader.gettingOnlyCAs()
        
        self.coordsets = reader.read()
        
        self.number_of_conformations = reader.numberOfFrames
        self.number_of_atoms = reader.numberOfAtoms
        
        self.pdb1_number_of_conformations = get_number_of_frames(self.pdb1)
        self.pdb2_number_of_conformations = get_number_of_frames(self.pdb2)
    
    def do_pdb_dic(self,pdb):
        return {
                  "source":pdb,
                  "conformations": get_number_of_frames(pdb)
        }
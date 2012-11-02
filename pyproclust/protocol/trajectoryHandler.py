'''
Created on 19/09/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
from pyproclust.tools.pdbTools import getPDBStructure, get_number_of_frames

class TrajectoryHandler(object):

    def __init__(self, pdb1, pdb2, rmsd_selection, subset = 'calpha'):
        if pdb1 == "" and pdb2 == "": 
            common.print_and_flush( "[ERROR] no pdb file(s) were given for RMSD calculation. Exiting...\n")
            exit()
        else:
            self.pdb1 = pdb1
            self.pdb2 = pdb2
            self.pdb_structure = getPDBStructure(pdb1, pdb2, rmsd_selection, subset)
            self.coordsets = self.pdb_structure.getCoordsets()
            self.number_of_conformations = len(self.coordsets)
            self.number_of_atoms = len(self.pdb_structure)
            self.pdb1_number_of_conformations = get_number_of_frames(self.pdb1)
            self.pdb2_number_of_conformations = get_number_of_frames(self.pdb2)
        
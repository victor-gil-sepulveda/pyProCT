"""
Created on 08/07/2014

@author: victor
"""
import numpy
from pyRMSD.condensedMatrix import CondensedMatrix
import pyproct.tools.mathTools as mathTools 

class DihedralRMSDBuilder(object):
    """
    """
    def __init__(self):
        """
        Constructor
        """
        pass

    
    @classmethod
    def build(cls, structure_data):
        print "Calculating dihedral RMSD matrix.  This may take some time ..."
        all_dihedrals = []
        num_conformations = structure_data.get_num_elements()
        for i in range(num_conformations):
            all_dihedrals.append(structure_data.get_dihedrals_for_conformation(i))

        data = []
        all_dihedrals = numpy.array(all_dihedrals)
        for i in range(num_conformations-1):
            dihedrals_i = all_dihedrals[i]
            for j in range(i+1, num_conformations):
                data.append(mathTools.angular_rmsd(dihedrals_i,all_dihedrals[j]))
        return CondensedMatrix(numpy.array(data))


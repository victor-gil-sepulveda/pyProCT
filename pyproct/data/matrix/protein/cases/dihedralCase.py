"""
Created on 08/07/2014

@author: victor
"""
import numpy
import math
from prody.measure.measure import calcPhi, calcPsi
from pyRMSD.condensedMatrix import CondensedMatrix

def rmsd(a,b):
    return math.sqrt(((a-b)**2).sum()/len(a))

class DihedralRMSDMatrixCalculator(object):
    """
    classdocs
    """
    def __init__(self):
        """
        Constructor
        """
        pass

    @classmethod
    def calculateDihedralsForCoordset(cls,structure, coordset):
        structure.setCoords(coordset)
        dihedral_angles = []
        for residue in structure.iterResidues():
            try:
                dihedral_angles.append(calcPhi(residue, radian=False, dist=None))
            except ValueError:
                dihedral_angles.append(0)

            try:
                dihedral_angles.append(calcPsi(residue, radian=False, dist=None))
            except ValueError:
                dihedral_angles.append(0)

        # 0 links with Nth residue and Nth with 0th. Those values are not needed anyway.
        return dihedral_angles[1:-1]

    @classmethod
    def build(cls, structure):
        print "Calculating dihedral RMSD matrix.  This may take some time ..."
        all_dihedrals = []
        coordsets =  structure.getCoordsets()
        for coordset in coordsets:
            dihedral_angles =  cls.calculateDihedralsForCoordset( structure, coordset)
            all_dihedrals.append(dihedral_angles)

        data = []
        all_dihedrals = numpy.array(all_dihedrals)
        for i in range(structure.numCoordsets()-1):
            dihedrals_i = all_dihedrals[i]
            for j in range(i+1, structure.numCoordsets()):
                data.append(rmsd(dihedrals_i,all_dihedrals[j]))
#                 print "1",dihedrals_i.tolist()
#                 print "2",all_dihedrals[j].tolist()
#                 print "RMSD", data[-1]
        return CondensedMatrix(numpy.array(data))
#         return CondensedMatrix(data)


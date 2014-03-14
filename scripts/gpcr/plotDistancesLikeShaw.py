'''
Created on 13/03/2014

@author: victor
'''


import prody
import os.path
from pyRMSD.RMSDCalculator import RMSDCalculator
import numpy
import math
import matplotlib
import pylab


alpha_atoms_selection = "name CA"

data = [
    {
        'dir':'2rh1',
        'selection':{
            'motif': "backbone and resnum 322:327",
            'residues': "ca resnum 131 272",
        }
    }
]

cwd = os.getcwd()
for datum in data:
    prot_name = datum['dir']
    print "========================\nWorking with %s\n========================"%(prot_name)
    # Look for the directory and enter it
    base_dir = os.path.join(cwd, prot_name)
    os.chdir(base_dir)

    pdb = prody.parsePDB("%s.pdb"%prot_name)
    ca_coordsets = pdb.select(alpha_atoms_selection)
    motif_coordsets = pdb.select(datum['selection']['motif'])

    calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
                 fittingCoordsets = ca_coordsets.getCoordsets(),
                 calculationCoordsets = motif_coordsets.getCoordsets())

    motif_rmsd = calculator.oneVsTheOthers(  conformation_number = 0,
                                             get_superposed_coordinates = False)


    arg131_leu272 = pdb.select(datum['selection']['residues'])

    calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
                 fittingCoordsets = ca_coordsets.getCoordsets(),
                 calculationCoordsets = arg131_leu272.getCoordsets())

    residues_rmsd_array, rearranged_coords, residues_rearranged_coords = calculator.oneVsTheOthers(
                                             conformation_number = 0,
                                             get_superposed_coordinates = True)

    residue_distances = []
    for conf in residues_rearranged_coords:
        arg131 = conf[0]
        leu272 = conf[1]
        r = conf[0] - conf[1]
        residue_distances.append(math.sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]))

    prody.writePDB("ca_atoms", ca_coordsets)
    prody.writePDB("motif", motif_coordsets)
    prody.writePDB("residues", arg131_leu272)

    print len(motif_rmsd), len(residue_distances)

    matplotlib.pyplot.scatter(residue_distances[1:], motif_rmsd)

    matplotlib.pyplot.show()




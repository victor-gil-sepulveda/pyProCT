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
#     {
#         'dir':'2rh1_2',
#         'selection':{
#             'motif': "backbone and resnum 322:327",
#             'residues': "ca resnum 131 272",
#         }
#     },
#     {
#         'dir':'2rh1_3',
#         'selection':{
#             'motif': "backbone and resnum 322:327",
#             'residues': "ca resnum 131 272",
#         }
#      },
    {
        'dir':'2rh1_4_Nano',
        'pdb_traj':'_2rh1_4_Nano.pdb',
        'native':'2RH1_Inac_Nano.pdb',
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

    native =  prody.parsePDB("%s"%datum['native'], chain='A',subset='backbone')
    pdb = prody.parsePDB("%s"%datum['pdb_traj'], chain='A',subset='backbone')

    print native.getCoordsets().shape
    print pdb.getCoordsets().shape
    native_ca_coordset = native.select(alpha_atoms_selection)
    native_motif_coordset = native.select(datum['selection']['motif'])
    ca_coordsets = pdb.select(alpha_atoms_selection)
    motif_coordsets = pdb.select(datum['selection']['motif'])

    print native_ca_coordset.getCoordsets().shape,ca_coordsets.getCoordsets().shape
    print native_motif_coordset.getCoordsets().shape,motif_coordsets.getCoordsets().shape
    exit()
#     calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
#                  fittingCoordsets = ca_coordsets.getCoordsets(),
#                  calculationCoordsets = motif_coordsets.getCoordsets())
    calculator = RMSDCalculator( calculatorType = "QCP_OMP_CALCULATOR",
                 fittingCoordsets = motif_coordsets.getCoordsets())

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




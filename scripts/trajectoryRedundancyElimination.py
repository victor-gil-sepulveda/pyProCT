'''
Created on 13/09/2012

@author: victor
'''

"""
This script calculates the rmsd matrix of a trajectory and generates a new trajectory
where all pair have a rmsd of no more than 'threshold'.
"""

import sys
import pyRMSD.RMSD
import prody.proteins
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.tools.pdbTools import extract_frames_from_trajectory_sequentially
import pickle 
                                                                                                                                                
if __name__ == '__main__':
    trajectory_pdb = sys.argv[1]
    threshold = float(sys.argv[2])
    final_trajectory_pdb = sys.argv[3]
    if sys.argv[4] == "load":
        rmsd_file_handler = open("rmsd.bin","r")
        matrix_data = pickle.load(rmsd_file_handler)
        rmsd_file_handler.close()
        matrix = CondensedMatrix(list(matrix_data))
    else:
        # Getting coordinates
        pdb = prody.proteins.parsePDB(trajectory_pdb)
        print "PDB parsed (",trajectory_pdb,")"
        coordsets = pdb.getCoordsets()
        number_of_conformations = len(coordsets)
        number_of_atoms = len(pdb)
        
        print "Calculating rmsd..."
        # Calculating rmsd
        rmsd = pyRMSD.RMSD.calculateRMSDCondensedMatrix(coordsets, "OMP_CALCULATOR")
        matrix = CondensedMatrix(rmsd)
        # Save it
        rmsd_file_handler = open("rmsd.bin","w")
        pickle.dump(matrix.get_data(),rmsd_file_handler)
        rmsd_file_handler.close()
    
    print "Filtering..."
    # Now get rid of redundancy!
    conformations_to_eliminate = []
    for i in range(matrix.row_length-1):
        for j in range(i+1,matrix.row_length):
            if matrix[i,j]<=threshold:
                conformations_to_eliminate.append(j)
    final_conformations = list(set(range(matrix.row_length))-set(conformations_to_eliminate))
    print "Erasing ",len(final_conformations)," structures"
    print "Writing..."
    traj_handler_in = open(trajectory_pdb,"r")
    traj_handler_out = open(final_trajectory_pdb,"w")
    extract_frames_from_trajectory_sequentially(traj_handler_in,matrix.row_length,traj_handler_out,final_conformations)
    traj_handler_in.close()
    traj_handler_out.close()
    
    
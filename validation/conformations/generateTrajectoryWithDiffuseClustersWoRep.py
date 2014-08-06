"""
Created on 10/09/2012

@author: victor
"""

import sys
import prody
import numpy
import numpy.random
from validation.conformations.generateTrajectoryWithNaturalClusters import preprocess_pdb,\
    print_matrix

MIN_SIMILARITY = 0.2
MAX_SIMILARITY = 0.8
MAX_TRIES = 20

if __name__ == '__main__':
    pdb, input_coordsets, cluster_frames, output = preprocess_pdb(sys.argv)

    # Generate new coordsets for each cluster. Any new frame for a given starting conformation
    # will be an interpolation between this conformation and all the other intial confs.
    for i in range(0,len(cluster_frames)-1):
        base_coordinates = input_coordsets[i]
        numatoms = len(base_coordinates)
        for j in range(i+1,len(cluster_frames)):
            directions = input_coordsets[j] - base_coordinates

            step = (MAX_SIMILARITY-MIN_SIMILARITY) / MAX_TRIES
            for k in range(MAX_TRIES):
                pdb.addCoordset(numpy.array([((MIN_SIMILARITY+(step*k))*directions)+base_coordinates]))

    print_matrix(pdb.select("name CA").getCoordsets(), output+"_big")
    prody.writePDB(output+".pdb", pdb)

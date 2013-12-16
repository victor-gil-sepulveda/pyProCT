'''
Created on 10/09/2012

@author: victor
'''

import sys
import prody
import numpy
import numpy.random
from validation.conformations.generateTrajectoryWithNaturalClusters import preprocess_pdb,\
    print_matrix

MIN_SIMILARITY = 0.05
MAX_SIMILARITY = 0.6
MAX_TRIES = 20

if __name__ == '__main__':
    pdb, input_coordsets, cluster_frames, output = preprocess_pdb(sys.argv)

    # Generate new coordsets for each cluster. Any new frame for a given starting conformation
    # will be an interpolation between this conformation and all the other intial confs.
    for i in range(len(cluster_frames)):
        base_coordinates = input_coordsets[i]
        numatoms = len(base_coordinates)
        for j in range(len(cluster_frames)):
            if i!= j:
                directions = input_coordsets[j] - base_coordinates

                displacement_magnitudes = numpy.random.exponential(scale = 1, size = 20)
                displacement_magnitudes /= numpy.max(displacement_magnitudes)
                displacement_magnitudes *= (MAX_SIMILARITY-MIN_SIMILARITY)
                displacement_magnitudes += MIN_SIMILARITY

                step = (MAX_SIMILARITY-MIN_SIMILARITY) / MAX_TRIES
                for k in range(MAX_TRIES):
                    pdb.addCoordset(numpy.array([((MIN_SIMILARITY+(step*k))*directions.T).T+base_coordinates]))

    print_matrix(pdb.select("name CA").getCoordsets(), output+"_big")
    prody.writePDB(output+".pdb", pdb)

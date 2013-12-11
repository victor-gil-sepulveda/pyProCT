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

MAX_SIMILARITY = 0.45

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
                for k in range(20):
                    displacement_magnitudes = numpy.random.uniform(low=0.0, high=MAX_SIMILARITY, size = numatoms)
                    pdb.addCoordset(numpy.array([ (displacement_magnitudes *directions.T).T+base_coordinates]))

    print_matrix(pdb.select("name CA").getCoordsets(), output+"_big")
    prody.writePDB(output+".pdb", pdb)


#                 Normalize translations to have at most a 1A displacement
#                 max_norm = numpy.max(numpy.sqrt(numpy.sum(directions*directions,1)))
#                 norm_dirs = directions/max_norm
#                 norm_dirs = (directions.T / numpy.sqrt(numpy.sum(directions*directions,1))).T
#                 for k in range(20):
#                     # New direction will be a random fraction of the distance with random sense
#                     # plus a totally random fluctuation
#                     fluctuation = numpy.random.uniform(low=-1.0, high=1.0, size = (numatoms,3))
#                     direction_modif = numpy.random.uniform(low=-1.0, high=1.0, size = (numatoms,3))
#                     translations = (direction_modif*(directions+fluctuation))
#                     norm_translations =(translations.T / numpy.sqrt(numpy.sum(translations*translations,1))).T
#                     displacement_magnitudes = numpy.random.randint(low=0, high=MAX_TRANSLATION, size = numatoms)
#                     pdb.addCoordset(numpy.array([ (displacement_magnitudes *norm_translations.T).T+base_coordinates]))

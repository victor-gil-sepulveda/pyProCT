'''
Created on 10/09/2012

@author: victor
'''

import sys
import prody
import numpy
import numpy.random
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix
import pyRMSD.RMSD
import matplotlib.pyplot as plt

if __name__ == '__main__':
    pdb_path = sys.argv[1]
    output = sys.argv[2]
    output_image_matrix = sys.argv[3]
    cluster_frames_str = sys.argv[4:len(sys.argv)]
    cluster_frames = []
    for c in cluster_frames_str:
        cluster_frames.append(int(c))
    
    print "Using frames ",cluster_frames
    
    pdb = prody.parsePDB(pdb_path)
    
    coordsets = numpy.array(pdb.getCoordsets())
    
    pdb.delCoordset(range(len(coordsets)))
    
    # Generate new coordsets for each cluster
    # Any new frame for a given starting conformation
    # will be an interpolation between this conformation and 
    # other of the chosen conformations.
    
    for i in cluster_frames:
        base_coordinates = coordsets[i]
        numatoms = len(base_coordinates)
        for j in cluster_frames:
            if i!= j:
                directions = coordsets[j] - base_coordinates
                for k in range(20):
                    # New direction will be a random fraction of the distance
                    # provoking some overlap
                    rand_list_array = numpy.random.exponential(5,numatoms*3)
                    # Normalize random numbers between [-0,6,0.6]
                    max_value = numpy.max(numpy.abs(rand_list_array))
                    norm_rand = rand_list_array  / max_value
                    dir_modif = numpy.resize(norm_rand,(numatoms,3))
                    pdb.addCoordset(numpy.array([dir_modif*directions+base_coordinates]))
    
    prody.writePDB(output, pdb)
    rmsd = pyRMSD.RMSD.calculateRMSDCondensedMatrix(pdb.getCoordsets(), "OMP_CALCULATOR")

    condensed_distance_matrix = CondensedDistanceMatrix(rmsd)
    
    # Normalize
    contents = condensed_distance_matrix.get_data()
    _max = numpy.max(contents)
    _min = numpy.min(contents)
    
    norm_contents = (contents - _min) / (_max - _min)
    del condensed_distance_matrix
    
    norm_condensed = CondensedDistanceMatrix(norm_contents)
    _max = numpy.max(norm_contents)
    _min = numpy.min(norm_contents)
    print _max," ",_min
    complete = numpy.zeros([norm_condensed.row_length,norm_condensed.row_length] ,dtype=numpy.float)
    
    for i in range(norm_condensed.row_length-1):
        for j in range(i+1,norm_condensed.row_length):
            complete[i][j] = norm_condensed[i,j]
            complete[j][i] = norm_condensed[i,j]
    
    fig = plt.figure()
    plt.gray()
    ax = fig.add_subplot(111)
    cax = ax.imshow(complete, interpolation='nearest')
    plt.show()
    plt.savefig(output_image_matrix)
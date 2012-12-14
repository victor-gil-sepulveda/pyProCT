'''
Created on 16/08/2012

@author: victor
'''
import sys
import prody
import random
import numpy

if __name__ == '__main__':
    pdb_path = sys.argv[1]
    output = sys.argv[2]
    cluster_frames_str = sys.argv[3:len(sys.argv)]
    cluster_frames = []
    for c in cluster_frames_str:
        cluster_frames.append(int(c))
    
    print "Using frames ",cluster_frames
    
    pdb = prody.parsePDB(pdb_path)
    
    coordsets = numpy.array(pdb.getCoordsets())
    
    pdb.delCoordset(range(len(coordsets)))
    
    # Generate new coordsets for each cluster
    for f in cluster_frames:
        base_coordinates = coordsets[f]
        numatoms = len(base_coordinates)
        for i in range(100):
            # Generate a random vector
            rand_list = []
            for r in range(numatoms*3):
                rand_list.append(random.uniform(-1, 1))
            coordset = numpy.resize(numpy.array(rand_list),(numatoms,3))
            pdb.addCoordset(numpy.array([coordset+base_coordinates]))
    prody.writePDB(output, pdb)
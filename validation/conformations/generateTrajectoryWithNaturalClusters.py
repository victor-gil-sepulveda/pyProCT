"""
Created on 16/08/2012

@author: victor
"""

import sys
import prody
import random
import numpy
from pyproct.tools.scriptTools import create_directory
from pyRMSD.RMSDCalculator import RMSDCalculator
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproct.tools.plotTools import matrixToImage

def get_frame_numbers(args):
    cluster_frames_str = sys.argv[3:len(args)]
    cluster_frames = []
    for c in cluster_frames_str:
        cluster_frames.append(int(c))
    cluster_frames = numpy.array(cluster_frames)
    print "Using frames ",cluster_frames
    return cluster_frames

def preprocess_pdb(args):
    pdb_file = args[1]
    output = "./" + args[2]+"/"+args[2]
    create_directory("./" + args[2])
    cluster_frames = get_frame_numbers(args)
    pdb = prody.parsePDB(pdb_file)
    # Get a copy of the pdb coords
    input_coordsets = numpy.array(pdb.getCoordsets()[cluster_frames])

    # Empty pdb
    pdb.delCoordset(range(pdb.numCoordsets()))

    # Build another pdb to store it
    input_pdb = prody.parsePDB(pdb_file)
    input_pdb.delCoordset(range(input_pdb.numCoordsets()))
    # And add the chosen coordsets
    for i in range(len(cluster_frames)):
        input_pdb.addCoordset(input_coordsets[i])
    prody.writePDB(output+"_ini.pdb", input_pdb)
    print_matrix(input_pdb.select("name CA").getCoordsets(), output)
    return pdb, input_coordsets, cluster_frames, output

def print_matrix(input_coordsets, output):
    # Generate the matrix and print it
    calculator = RMSDCalculator(calculatorType="QCP_OMP_CALCULATOR", fittingCoordsets = input_coordsets)
    matrixToImage(CondensedMatrix(calculator.pairwiseRMSDMatrix()), output + ".png")

if __name__ == '__main__':
    pdb, input_coordsets, cluster_frames, output = preprocess_pdb(sys.argv)

    # Generate new coordsets for each cluster
    for f in range(len(cluster_frames)):
        base_coordinates = input_coordsets[f]
        numatoms = len(base_coordinates)
        for i in range(100):
            # Generate a random vector
            rand_list = []
            for r in range(numatoms*3):
                rand_list.append(random.uniform(-1, 1))
            rand_displacements = numpy.resize(numpy.array(rand_list),(numatoms,3))
            pdb.addCoordset(numpy.array([rand_displacements+base_coordinates]))
    print_matrix(pdb.select("name CA").getCoordsets(), output+"_big")
    prody.writePDB(output+".pdb", pdb)

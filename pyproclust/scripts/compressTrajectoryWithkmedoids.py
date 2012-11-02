#!/usr/bin/python

'''
Created on 14/09/2012

@author: victor
'''

import sys
import optparse
import prody.proteins
from pyproclust.protocol.matrixHandler import MatrixHandler
from pyproclust.algorithms.kmedoids.kMedoids import KMedoids
from pyproclust.tools.pdbTools import extract_frames_from_trajectory
import random
import math

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog -m <arg> -c <arglist> [-o <arg>]', version='1.0')
    
    parser.add_option('-k', action="store", type='int', dest = "k", help="Number of clusters to get.",  metavar = "1")
    parser.add_option('-o', action="store", dest = "pdb_out", help="Path and name for the output trajectory.",metavar = "out.pdb")
    parser.add_option('-i', action="store", dest = "pdb_in", help="Path and name for the input trajectory.",metavar = "trajectory.pdb")
    parser.add_option("-w", action="store_true", dest="population_weighted")
    parser.add_option("-r", action="store_true", dest="random_centroid")
    parser.add_option('-n', action="store", type='int', dest = "n", help="Number of structures to get (wen using -w -r).",  metavar = "50")
    options, args = parser.parse_args()
    
    if not (options.k):
        parser.error("Please specify the number of clusters.")
    
    if not (options.pdb_out):
        parser.error("Please specify the output trajectory file.")
    
    if not (options.pdb_in):
        parser.error("Please specify the input trajectory file.")
        
    print "Parsing PDB : ",options.pdb_in 
    sys.stdout.flush()
    pdb = prody.proteins.parsePDB(options.pdb_in)
    coordsets = pdb.getCoordsets()[100:]
    number_of_conformations = len(coordsets)
    number_of_atoms = len(pdb)
    selection = pdb.select("backbone or not protein")
    pdb = selection.copy()
    
    print "Calculating RMSD matrix ..."
    sys.stdout.flush()
    matrixHandler = MatrixHandler()
    matrixHandler.createMatrix(coordsets)
    matrixHandler.saveMatrix("matrix")
    
    print "Clustering ..." 
    sys.stdout.flush()
    kmedClustering = KMedoids(matrixHandler.distance_matrix)
    clustering = kmedClustering.perform_clustering({"k":options.k,"seeding_max_cutoff":8.})
    
    print "Getting centers ..."
    sys.stdout.flush()
    centers = []
    populations = []
    if options.random_centroid == False:
        for cluster in clustering.clusters:
            centers.append(cluster.calculate_medoid(matrixHandler.distance_matrix))
            populations.append(cluster.get_size())
    else:
        if options.population_weighted == True:
            # Getting number of structures we want of each of the clusters and this structures
            for cluster in clustering.clusters:
                if len(centers)< options.n:# trying to avoid rounding errors
                    num_structs = int(math.ceil((cluster.get_size() / float(number_of_conformations))*options.n))
                    if len(centers) + num_structs >= options.n:
                        num_structs = options.n - len(centers)
                    print "We'll get ", num_structs, "random conformations from cluster (",cluster.get_size(),")"
                    temporary_list = list(cluster.all_elements)
                    random.shuffle(temporary_list)
                    centers.extend(temporary_list[0:num_structs])
                    populations.append(cluster.get_size())
                    print "Which are", len(centers), centers
    
    print "Saving..."
    sys.stdout.flush()
    traj_handler_in = open(options.pdb_in,"r")
    traj_handler_out = open(options.pdb_out,"w")
    extract_frames_from_trajectory(traj_handler_in,matrixHandler.distance_matrix.row_length,traj_handler_out, centers, "MODEL","ENDMDL",keep_header = True)
    traj_handler_in.close()
    traj_handler_out.close()
    
    print "Saving populations"
    pop_handler = open("pops","w")
    for p in populations:
        pop_handler.write(str(p)+"\n")
    pop_handler.close()
    
    
    
    
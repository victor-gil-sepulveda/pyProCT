#!/usr/local/bin/python

'''
Created on 09/02/2012

@author: victor
'''
import optparse
import pickle
import pyproclust.tools.commonTools as common
from pyproclust.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
from pyproclust.algorithms.hierarchical.hierarchicalTools import calculate_cutoff_numcluster_list

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog ', version='1.0')
    
    parser.add_option('--max-distance', action="store",type = 'float', dest = "max_dist",help="Maximum distance value in the distance matrix (upper bound)",metavar = "4.5")
    parser.add_option('-l', action="store", dest = "stored_linkage_matrix_file",help="A previously stored linkage matrix in binary form.",metavar = "link_matrix.bin")
    parser.add_option('-o', action="store", dest = "output_file",help="output file containing the detailed cutoff-cluster_number pairings and the plot (without extension)",metavar = "hierarchical_analysis")
        

    options, args = parser.parse_args()
    
    
    # Mandatory options
    if not (options.max_dist or options.stored_linkage_matrix_file):
        parser.error("Please see usage.")

   
    ########################################
    # Script Start
    ######################################## 
    common.print_and_flush("Loading linkage matrix...")  
    link_matrix_file_handler = open(options.stored_linkage_matrix_file,"r")
    linkage_matrix = pickle.load(link_matrix_file_handler)
    link_matrix_file_handler.close()
    common.print_and_flush(" Done\n")
    
    algorithm = HierarchicalClusteringAlgorithm(None)
    
    ##########################################
    # Find the threshold where the algorithm 
    # starts to retrieve clusters.
    ########################################## 
    common.print_and_flush("Calculating cutoffs ... ")  
    cutoffs = calculate_cutoff_numcluster_list(algorithm, linkage_matrix, options.max_dist)
    common.print_and_flush(" Done\n")

    ###########################
    # Save and show plot
    ###########################
    if options.output_file:
        cutoffs.sort()
        # Print it into a file
        output_handler = open(options.output_file+'_values.txt','w')
        for c in cutoffs:        
            output_handler.write(str(c[0])+" "+str(c[1])+"\n")
        output_handler.close()

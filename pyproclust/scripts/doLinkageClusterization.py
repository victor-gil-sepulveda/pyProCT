#!/usr/local/bin/python

'''
Created on 09/02/2012

@author: victor
'''
import optparse
import pickle
import pyproclust.tools.commonTools as common
import sys
from pyproclust.tools import scriptTools
from pyproclust.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator


def evaluate_clustering_for_all_cutoffs(algorithm, options,analyzer):
    
    for cutoff in options.cutofflist:
        common.print_and_flush( "Hierarchical clustering for cutoff "+str(cutoff)+" ...") 
        clusterization = algorithm.perform_clustering({"cutoff":cutoff})
        common.print_and_flush(" Done\n")
        if options.cluster_save_path:
            scriptTools.save_clusters_as_binary(options.cluster_save_path+"_cutoff_"+str(cutoff)+".bin",clusterization)
    
        common.print_and_flush("Running analysis ...")
        analyzer.run_analysis_for(clusterization)
        common.print_and_flush(" Done\n")
    
    
    if options.linkage_matrix_file:
        common.print_and_flush("Saving linkage matrix ...")
        file_handler = open(options.linkage_matrix_file+".bin","w")
        pickle.dump(algorithm.hie_mat, file_handler)
        file_handler.close()
        common.print_and_flush(" Done\n")

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog -m|-l <arg> -c <arglist> [-s <arg>] [-o <arg>]', version='1.0')
    
    parser.add_option('-m', action="store", dest = "matrix_file",help="The upper triangular matrix of the distance matrix in text form.",metavar = "matrix_rmsd.txt")
    parser.add_option('-l', action="store", dest = "stored_linkage_matrix_file",help="A previously stored linkage matrix.",metavar = "link_matrix.bin")
    parser.add_option('-s', action="store", dest = "linkage_matrix_file", help="Path of the file to store the linkage matrix (without extension).",metavar = "link_matrix")
    parser.add_option('-c', '--cutoffs', action="callback",dest = "cutofflist", callback =  common.vararg_callback ,metavar = "1. 2. 3. ...")
    parser.add_option('-o', action="store", dest = "output", help="Path and name for the output file.",metavar = "out.txt")
    parser.add_option('--clusters', action="store",dest = "cluster_save_path", help="Prefix for the clustering binary file which will be stored inside ./clusterings",metavar = "clusters.bin")
    
    options, args = parser.parse_args()
    
    # Mandatory options
    if not (options.matrix_file or options.stored_linkage_matrix_file):
        parser.error("Please specify the matrix file.")
    if len(options.cutofflist) == 0:
        parser.error("Please specify the cutoff list.")

   
    ########################################
    # Script Start
    ######################################## 
    linkage_matrix = None
    if options.matrix_file:
        ########################################
        # Load distance matrix
        ######################################## 
        common.print_and_flush("Loading condensed matrix...")
        condensed_distance_matrix = scriptTools.load_matrix(options.matrix_file)
        common.print_and_flush(" Done\n")
    
    if options.stored_linkage_matrix_file:
        ########################################
        # Load linkage matrix from disk
        ########################################
        common.print_and_flush("Loading linkage matrix...")  
        link_matrix_file_handler = open(options.stored_linkage_matrix_file,"r")
        linkage_matrix = pickle.load(link_matrix_file_handler)
        link_matrix_file_handler.close()
        common.print_and_flush(" Done\n")
    
    ############################################################
    # Create the hierarchical algorithm and perform evaluation
    ############################################################
    algorithm = HierarchicalClusteringAlgorithm(condensed_distance_matrix)
    algorithm.hie_mat = linkage_matrix
    
    ##########################################
    # Evaluate the algorithm for all the cutoffs
    ########################################## 
    analyzer = AnalysisRunner()
    AnalysisPopulator(analyzer, condensed_distance_matrix)
    evaluate_clustering_for_all_cutoffs(algorithm,options,analyzer)
    
    ########################################
    # Report writing
    ########################################
    if options.output:
        scriptTools.save_report(options.output,analyzer)
    else:
        sys.stdout.write(analyzer.generate_report())
    
    

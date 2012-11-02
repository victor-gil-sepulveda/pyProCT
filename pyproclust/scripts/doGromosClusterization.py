#!/usr/local/bin/python
'''
Created on 28/02/2012

@author: victor
'''
import optparse
import pyproclust.tools.commonTools as common
import sys
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator
from pyproclust.tools import scriptTools
from pyproclust.algorithms.gromos.gromosAlgorithm import GromosAlgorithm
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner

def evaluate_clustering_for_all_cutoffs( algorithm, options,analyzer):
    for cutoff in options.cutofflist:
        common.print_and_flush( "GROMOS clustering for cutoff "+str(cutoff)+" ...") 
        clusterization = algorithm.perform_clustering({"cutoff":cutoff})
        if options.del_noise:
            clusterization.eliminate_noise(options.del_noise)
        if options.cluster_save_path:
            scriptTools.save_clusters_as_binary(options.cluster_save_path+"_cutoff_"+str(cutoff)+".bin",clusterization)
        common.print_and_flush(" Done\n")
        common.print_and_flush("Running analysis ...")
        analyzer.run_analysis_for(clusterization)
        common.print_and_flush(" Done\n")

if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage='%prog -m <arg> -c <arglist> [-o <arg>]', version='1.0')
    
    parser.add_option('-m', action="store", dest = "matrix_file",help="The upper triangular matrix of the distance matrix in text form.",metavar = "matrix_rmsd.txt")
    parser.add_option('-c', '--cuttoffs', action="callback",dest = "cutofflist", callback =  common.vararg_callback ,metavar = "1. 2. 3. ...")
    parser.add_option('-o', action="store", dest = "output", help="Path and name for the output file.",metavar = "out.txt")
    parser.add_option('--del-noise', action="store", type='int',dest = "del_noise", help="Deletes clusters which size is lower than this.",metavar = "4")
    parser.add_option('--clusters', action="store",dest = "cluster_save_path", help="Prefix for the clustering binary file which will be stored inside ./clusterings",metavar = "gromos_clusters")
    
    options, args = parser.parse_args()
    
    
    # Mandatory options
    if not (options.matrix_file):
        parser.error("Please specify the matrix file.")
        
    if len(options.cutofflist) == 0:
        parser.error("Please specify one or more cutoffs (example '-c 1. 2. 3.' )")
         
    # More possible options: 
    # - Generate node graph (free energy, transition probability)
    # - Get information from the graph
    
    ########################################
    # Script Start
    ######################################## 
    
    
    ########################################
    # Load distance matrix
    ######################################## 
    common.print_and_flush("Loading condensed matrix...")
    condensed_distance_matrix = scriptTools.load_matrix(options.matrix_file)
    common.print_and_flush(" Done\n")
    
    ########################################
    # Algorithm creation
    ########################################
    algorithm = GromosAlgorithm(condensed_distance_matrix)
    
    ########################################
    # Clusterization / Evaluation
    ########################################
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

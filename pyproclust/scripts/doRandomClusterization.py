#!/usr/local/bin/python
'''
Created on 28/02/2012

@author: victor
'''
import optparse
import pyproclust.tools.commonTools as common
import sys
from pyproclust.tools import scriptTools
from pyproclust.algorithms.random.RandomAlgorithm import RandomClusteringAlgorithm
from pyproclust.algorithms.random.FakeDistributionRandomAlgorithm import FakeDistributionRandomClusteringAlgorithm
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator

def evaluate_clustering_for_numclust_or_distrib( algorithm, options,analyzer):
    
    if options.max_clusters != None:
        for max_clusters in options.max_clusters:
            common.print_and_flush( "Random clustering ...")
            clusterization = algorithm.perform_clustering({"max_num_of_clusters":max_clusters})
            common.print_and_flush(" Done\n")
            if options.cluster_save_path:
                scriptTools.save_clusters_as_binary(options.cluster_save_path+"_"+str(max_clusters),clusterization)
            common.print_and_flush("Running analysis ...")
            analyzer.run_analysis_for(clusterization)
            common.print_and_flush(" Done\n")
    else:
        common.print_and_flush( "Random distribution clustering ...")
        clusterization = algorithm.perform_clustering({"distribution":options.rand_distr})
        common.print_and_flush(" Done\n")
        if options.cluster_save_path:
            scriptTools.save_clusters_as_binary(options.cluster_save_path,clusterization)
        common.print_and_flush("Running analysis ...")
        analyzer.run_analysis_for(clusterization)
        common.print_and_flush(" Done\n")

if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage='%prog -m <arg> --maxclusters <arglist> [-o <arg>]', version='1.0')
    
    parser.add_option('-m', action="store", dest = "matrix_file",help="The upper triangular matrix of the distance matrix in text form.",metavar = "matrix_rmsd.txt")
    parser.add_option('--maxclusters', action="callback",dest = "max_clusters", callback =  common.vararg_callback ,metavar = "1 2 3 ...")
    parser.add_option('--random-distribution', action="callback",dest = "rand_distr", callback =  common.vararg_callback ,metavar = "60 33 7")
    parser.add_option('-o', action="store", dest = "output", help="Path and name for the output file.",metavar = "out.txt")
    parser.add_option('--clusters', action="store",dest = "cluster_save_path", help="Prefix for the clustering binary file which will be stored inside ./clusterings",metavar = "clusters.bin")
    
    options, args = parser.parse_args()
    
    
    # Mandatory options
    if not (options.matrix_file or options.stored_linkage_matrix_file):
        parser.error("Please specify the matrix file.")
        
    if options.max_clusters== None and options.rand_distr == None:
        parser.error("Either maxclusters or random-distribution are mandatory.")
    
    # Mutual exclusion
    if options.max_clusters != None and options.rand_distr != None:
        parser.error("Cannot use both modes at the same time.")
    
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
    if options.max_clusters != None:
        algorithm = RandomClusteringAlgorithm(condensed_distance_matrix)
    else:
        algorithm = FakeDistributionRandomClusteringAlgorithm(condensed_distance_matrix)
    
    ########################################
    # Clusterization / Evaluation
    ########################################
    analyzer = AnalysisRunner()
    AnalysisPopulator(analyzer, condensed_distance_matrix)
    evaluate_clustering_for_numclust_or_distrib(algorithm,options,analyzer)
    
    ########################################
    # Report writing
    ########################################
    if options.output:
        scriptTools.save_report(options.output,analyzer)
    else:
        sys.stdout.write(analyzer.generate_report())
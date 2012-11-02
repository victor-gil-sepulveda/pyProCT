#!/usr/local/bin/python
'''
Created on 18/04/2012

@author: victor
'''
import optparse
import pyproclust.tools.commonTools as common
import sys
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator
from pyproclust.tools import scriptTools
from pyproclust.algorithms.dbscan.dbscanAlgorithm import DBSCANAlgorithm
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner


def evaluate_clustering(algorithm,options,analyzer):
    common.print_and_flush( "DBSCAN clustering for eps "+str(options.eps)+" and minpts "+str(options.minpts)+" ...") 
    clusterization = algorithm.perform_clustering({"eps":options.eps,"minpts":options.minpts})
    common.print_and_flush(" Done\n")
    if options.cluster_save_path:
            scriptTools.save_clusters_as_binary(options.cluster_save_path+"_eps_"+str(options.eps)+"_minpts_"+str(options.minpts),clusterization)
    common.print_and_flush("Running analysis ...")
    analyzer.run_analysis_for(clusterization)
    common.print_and_flush(" Done\n")

if __name__ == '__main__':
    ## TODO change usage
    parser = optparse.OptionParser(usage='%prog -m <arg> --eps <arg> --minpts <arg> [-o <arg>] ',
                                   version='1.0')
    
    parser.add_option('-m', action="store", dest = "matrix_file", help="The upper triangular matrix of the distance matrix in text form.",metavar = "matrix_rmsd.txt")
    parser.add_option('--eps', type="float",action="store", dest = "eps", help="Eps, cutoff distance to get the neighbours of an element.",metavar = "2.0")
    parser.add_option('--minpts', type="int",action="store", dest = "minpts", help="Number of neighbors for a point to be considered a core point.",metavar = "4")
    parser.add_option('-o', action="store", dest = "output", help="Path and name for the output file.",metavar = "out.txt")
    parser.add_option('--clusters', action="store",dest = "cluster_save_path", help="Prefix for the clustering binary file which will be stored inside ./clusterings",metavar = "clusters.bin")
    
    options, args = parser.parse_args()
    
    # Mandatory options
    if not (options.matrix_file ):
        parser.error("Please specify the matrix file.")
    if not (options.eps):
        parser.error("Please specify the eps value.")
    if not (options.minpts):
        parser.error("Please specify the numpts value.")
        
        
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
    algorithm = DBSCANAlgorithm(condensed_distance_matrix)
    
    ########################################
    # Clusterization / Evaluation
    ########################################
    analyzer = AnalysisRunner()
    AnalysisPopulator(analyzer, condensed_distance_matrix)
    evaluate_clustering(algorithm,options,analyzer)
    
    ########################################
    # Report writing
    ########################################
    if options.output:
        scriptTools.save_report(options.output,analyzer)    
    else:
        sys.stdout.write(analyzer.generate_report())
        

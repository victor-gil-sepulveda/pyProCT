#!/usr/local/bin/python
'''
Created on 28/02/2012

@author: victor
'''
import optparse
import pyproclust.tools.commonTools as common
import sys
from pyproclust.clustering.analysis.analysisPopulator import AnalysisPopulator
from pyproclust.algorithms.kmedoids.kMedoids import KMedoids
from pyproclust.tools import scriptTools
from pyproclust.clustering.analysis.analysisRunner import AnalysisRunner

def evaluate_clustering_for_all_cutoffs( algorithm, options,analyzer):
    for k in options.klist:
        common.print_and_flush( "Kmedoids clustering for k "+str(k)+" ...") 
        clusterization = algorithm.perform_clustering({"k":k,"seeding_max_cutoff":2.0})
        common.print_and_flush(" Done\n")
        common.print_and_flush("Running analysis ...")
        analyzer.run_analysis_for(clusterization)
        common.print_and_flush(" Done\n")

if __name__ == '__main__':
    
    parser = optparse.OptionParser(usage='%prog -m <arg> -c <arglist> [-o <arg>]', version='1.0')
    
    parser.add_option('-m', action="store", dest = "matrix_file",help="The upper triangular matrix of the distance matrix in text form.",metavar = "matrix_rmsd.txt")
    parser.add_option('-k',  action="callback",dest = "klist", callback =  common.vararg_callback ,metavar = "1. 2. 3. ...")
    parser.add_option('-o', action="store", dest = "output", help="Path and name for the output file.",metavar = "out.txt")
    
    options, args = parser.parse_args()
    
    
    # Mandatory options
    if not (options.matrix_file):
        parser.error("Please specify the matrix file.")
        
    if len(options.klist) == 0:
        parser.error("Please specify one or more cutoffs (example '-k 1 2 3' )")
         
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
    algorithm = KMedoids(condensed_distance_matrix)
    
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

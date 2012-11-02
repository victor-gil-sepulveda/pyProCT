#!/usr/local/bin/python
'''
Created on 03/05/2012

@author: victor
'''
import pyproclust.tools.commonTools as common
import optparse
import scipy.stats
import Gnuplot
import numpy
import os
from pyproclust.algorithms.dbscan.dbscanTools import kdist
from pyproclust.tools import scriptTools
from pyproclust.tools.distanceMatrixAnalysisTools import kdist_analysis

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog -m <arg> [-dd <arg> -mdd <arg> --mean-dist-cutoff <arg> --elem-percent-cutoff <arg> --kdist <arg> ] ',
                                   version='1.0')
    
    parser.add_option('-m', action="store", dest = "matrix_file", help="The upper triangular matrix of the distance matrix in text form.",metavar = "matrix_rmsd.txt")
    parser.add_option('--dd', action="store_true",dest = "dist_ditrib", )
    parser.add_option('--mdd', action="store_true",dest = "mean_dist_distrib" )
    parser.add_option('--mean-dist-cutoff', action="callback",dest = "mean_dist", callback =  common.vararg_callback ,metavar = "1. 2. 3. ...")
    parser.add_option('--elem-percent-cutoff', action="callback",dest = "elem_percent", callback =  common.vararg_callback ,metavar = "1. 2. 3. ...")
    parser.add_option('--kdist', action="callback",dest = "kdistlist", callback =  common.vararg_callback ,metavar = "first last increment")
    
    options, args = parser.parse_args()
    
    # Mandatory options
    if not (options.matrix_file ):
        parser.error("Please specify the matrix file.")
    
    ########################################
    # Script Start
    ######################################## 
    
    ########################################
    # Load distance matrix
    ########################################
    common.print_and_flush("Loading condensed matrix...")
    condensed_distance_matrix=scriptTools.load_matrix(options.matrix_file)
    common.print_and_flush(" Done\n")
    
    ########################################
    # Do some analysis
    ########################################
    minimum, maximum = condensed_distance_matrix.get_minimum_and_maximum()

    print "--------------"
    print "Matrix values" 
    print "-------------"
    print "Minimum:\t","%.5f"%minimum
    print "Maximum:\t","%.5f"%maximum
    print "Mean:\t","%.5f"%(numpy.mean(condensed_distance_matrix.get_data()))
    print "Std. Dev.:\t","%.5f"%(numpy.std(condensed_distance_matrix.get_data()))
    print "Skewness:\t", "%.5f"%(scipy.stats.skew(condensed_distance_matrix.get_data()))
    print "Kurtosis:\t","%.5f"%(scipy.stats.kurtosis(condensed_distance_matrix.get_data()))

    ## Recreate matrix_images folder
    scriptTools.make_directory('./matrix_images')
    
    print "--------------"
    print "Matrix plots" 
    print "-------------"
    
    g = Gnuplot.Gnuplot(debug=1)
    # Distance distribution
    if options.dist_ditrib :
        g('set data style boxes')
        g.title('Distance distribution')
        g.plot(condensed_distance_matrix.distance_distribution(granularity = 50))
        g.hardcopy('./matrix_images/distance_distrib.ps', enhanced=1, color=1)
    #    raw_input('Please press return to continue...\n')
        g.reset()

    # Mean distance distribution
    if options.mean_dist_distrib :
        g('set data style lines')
        g.title('Mean distance distribution')
        g.plot(condensed_distance_matrix.mean_distance_distribution())
        g.hardcopy('./matrix_images/mean_distance_distrib.ps', enhanced=1, color=1)
    #    raw_input('Please press return to continue...\n')
        g.reset()
    
    
    # Element neighbors for cutoff
    if options.elem_percent:
        for c in options.elem_percent:
            g('set data style lines')
            g.title(' % of number of elements in cutoff = '+str(c))
            g.plot(condensed_distance_matrix.percent_of_elements_within_cutoff_per_element(c))
            g.hardcopy('./matrix_images/element_percent_per_element_'+str(c)+'.ps', enhanced=1, color=1)
#           raw_input('Please press return to continue...\n')
            g.reset()
        
    # Element mean distance for cutoff
    if options.mean_dist:
        for c in options.mean_dist:
            g('set data style lines')
            g.title('Mean distance plot for cutoff = '+str(c))
            g.plot(condensed_distance_matrix.mean_distance_per_element(c))
            g.hardcopy('./matrix_images/mean_dist_per_element_'+str(c)+'.ps', enhanced=1, color=1)
            g.reset()
        
    # K-dist
    if options.kdistlist:
        kdist_analysis(condensed_distance_matrix, options.kdistlist, './matrix_images/k-dist_')
        
    # Making the actual images (transforming ps files)
    os.system("find ./matrix_images -name '*.ps' -exec convert -background white -flatten -rotate 90 -resample 150 {} {}.png \;")
    
    if options.elem_percent:
        scriptTools.tile_images_using_montage('./matrix_images/element_percent_per_element_', 2, len(options.elem_percent),"./matrix_images/tiled_element_percent_per_element.png")
    
    if options.mean_dist:
        scriptTools.tile_images_using_montage('./matrix_images/mean_dist_per_element_', 2, len(options.mean_dist),"./matrix_images/tiled_mean_dist_per_element.png")
    
    if options.kdistlist:
        scriptTools.tile_images_using_montage('./matrix_images/k-dist_', 2, len(options.kdistlist),"./matrix_images/tiled_kdist.png")
    
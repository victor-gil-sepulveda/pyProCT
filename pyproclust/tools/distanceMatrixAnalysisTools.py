'''
Created on 30/05/2012

@author: victor
'''

import Gnuplot
from pyproclust.algorithms.dbscan.dbscanTools import kdist

def distance_distribution_analysis(condensed_distance_matrix,folder_path):
    g = Gnuplot.Gnuplot()
    # Distance distribution
    g('set data style boxes')
    g.title('Distance distribution')
    g.plot(condensed_distance_matrix.distance_distribution(granularity = 50))
    g.hardcopy(folder_path+'/distance_distrib.ps', enhanced=1, color=1)
    g.reset()

def mean_distance_distribution_analysis(condensed_distance_matrix,folder_path):
    g = Gnuplot.Gnuplot()
    g('set data style lines')
    g.title('Mean distance distribution')
    g.plot(condensed_distance_matrix.mean_distance_distribution())
    g.hardcopy(folder_path+'/mean_distance_distrib.ps', enhanced=1, color=1)
    g.reset()

def percent_of_elements_within_cutoff_analysis(condensed_distance_matrix,folder_path,cutoff):
    g = Gnuplot.Gnuplot()
    g('set data style lines')
    g.title(' % of number of elements in cutoff = '+str(cutoff))
    g.plot(condensed_distance_matrix.percent_of_elements_within_cutoff_per_element(cutoff))
    g.hardcopy(folder_path+'/element_percent_per_element_'+str(cutoff)+'.ps', enhanced=1, color=1)
    g.reset()

def mean_distance_per_element(condensed_distance_matrix,folder_path,cutoff):
    g = Gnuplot.Gnuplot()
    g('set data style lines')
    g.title('Mean distance plot for cutoff = '+str(cutoff))
    g.plot(condensed_distance_matrix.mean_distance_per_element(cutoff))
    g.hardcopy(folder_path+'/mean_dist_per_element_'+str(cutoff)+'.ps', enhanced=1, color=1)
    g.reset()

def kdist_analysis(condensed_distance_matrix,klist,folder_path = None):
    g = Gnuplot.Gnuplot()
    k_matrix = kdist(klist,condensed_distance_matrix)
    for i in len(klist):
        k = klist[i]
        k_array = k_matrix[i]
        g = Gnuplot.Gnuplot(debug=1)
        g('set data style lines')
        g.title('Kdist graph for k='+str(int(k)))
        g.plot(k_array)
        if folder_path:
            g.hardcopy(folder_path+str(int(k))+'.ps', enhanced=1, color=1)
        else:
            raw_input('Please press return to continue...\n')
        g.reset()

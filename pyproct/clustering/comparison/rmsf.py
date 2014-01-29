'''
Created on 19/11/2013

@author: victor
'''
import numpy
import math

def calc_rmsf(coordsets, reference):
    # Calculate distances, reference is prototype
    conf_minus_ref = coordsets-reference # Xi-Xref
    distances = numpy.sqrt(numpy.sum(conf_minus_ref**2,2))
    return math.sqrt(numpy.mean(distances**2, 0))

def calc_rmsf_of_cluster(alpha_carbons_trajectory_coordsets, cluster):
    if len(alpha_carbons_trajectory_coordsets) != 0:
        reference = alpha_carbons_trajectory_coordsets[cluster.prototype]
        cluster_coordsets = alpha_carbons_trajectory_coordsets[cluster.all_elements]
        return calc_rmsf(cluster_coordsets, reference)
    else:
        print "[WARNING][calc_rmsf_of_cluster]  No CA atoms found. Aborting operation."
        return

'''
Created on 19/11/2013

@author: victor
'''
import numpy

def calc_norms(coordsets, reference):
    # Calculate distances, reference is prototype
    conf_minus_ref = coordsets-reference
    norms = numpy.sqrt(numpy.sum(conf_minus_ref**2,2))
    return numpy.mean(norms,0)

def CA_mean_square_displacement_of_cluster(alpha_carbons_trajectory_coordsets, cluster):
    if len(alpha_carbons_trajectory_coordsets) != 0:
        reference = alpha_carbons_trajectory_coordsets[cluster.prototype]
        cluster_coordsets = alpha_carbons_trajectory_coordsets[cluster.all_elements]
        return calc_norms(cluster_coordsets, reference)
    else:
        print "[WARNING][CA_mean_square_displacement_of_cluster]  No CA atoms found. Aborting operation."
        return
    
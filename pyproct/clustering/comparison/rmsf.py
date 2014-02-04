'''
Created on 19/11/2013

@author: victor
'''
import numpy

def calc_rmsf_of_cluster(alpha_carbons_trajectory_coordsets, cluster):
    """
    Adapted from Prody code
    """
    if len(alpha_carbons_trajectory_coordsets) != 0:
        cluster_coordsets = alpha_carbons_trajectory_coordsets[cluster.all_elements]
        mean_conformation = cluster_coordsets.mean(0)
        ssqf = numpy.zeros(mean_conformation.shape)
        for conf in cluster_coordsets:
            ssqf += (conf - mean_conformation) ** 2
        return (ssqf.sum(1) / len(cluster_coordsets))**0.5
    else:
        print "[WARNING][calc_rmsf_of_cluster]  No CA atoms found. Aborting operation."
        return

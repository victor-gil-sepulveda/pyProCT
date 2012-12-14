'''
Created on 20/02/2012

@author: victor
'''

import scipy.spatial.distance as distance
import pyproclust.algorithms.validation.validationTools as vt
import pyproclust.algorithms.gromos.gromosAlgorithm as ga
import pyproclust.algorithms.dbscan.dbscanAlgorithm as dba
import pyproclust.algorithms.validation.datasets as data
from pyproclust.matrix.condensedMatrix import CondensedDistanceMatrix

"""
Script for visual validation of algorithms.
"""
if __name__ == '__main__':
    for dataset in data.all_datasets[0:4]:
        
        # Creating matrixes for datasets...
        observations = vt.dataset_loading_2D(dataset)
        total_number_of_observations = len(observations)
        condensed_matrix_data = distance.pdist(observations)
        condensed_matrix = CondensedDistanceMatrix(condensed_matrix_data)
        
        # Testing it for different algorithms...
        #####################
        ### GROMOS
        #####################
        algorithm = ga.GromosAlgorithm(condensed_matrix)
        clusterization =  algorithm.perform_clustering(4.0)
        clusterization.eliminate_noise(2)
        candidates = []
        
        for c in clusterization.clusters:
            if c.get_size()>3 and len(candidates)<5:
                candidates.append(c.prototype)
            
        vt.show_2D_dataset(observations,15 ,20 ,4. ,candidates).show()
        vt.show_2D_dataset_clusters(observations, 15, clusterization,20).show()
        #####################
        ### DBSCAN
        #####################
        algorithm = dba.DBSCANAlgorithm(condensed_matrix)
        clusterization =  algorithm.perform_clustering(4.0,3)
        vt.show_2D_dataset_clusters(observations, 15, clusterization,20).show()
        
        print "Done"
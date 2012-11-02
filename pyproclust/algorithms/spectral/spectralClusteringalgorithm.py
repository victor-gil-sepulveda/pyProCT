'''
Created on 14/08/2012

@author: victor
'''
from sklearn.cluster.spectral import SpectralClustering
from pyproclust.clustering.cluster import gen_clusters_from_class_list
from pyproclust.clustering.clusterization import Clustering

class SpectralClusteringAlgorithm(object):
    '''
    '''
    def __init__(self,condensed_matrix):
        self.condensed_matrix = condensed_matrix
        self.affinity_matrix = self.condensed_matrix.calculate_affinity_matrix(1)
    
    def perform_clustering(self,kwargs):
        # Mandatory parameter
        k = int(kwargs["k"])
#        {'k': 3, 'random_state': None, 'mode': 'amg', 'n_init': 10}
#        laplacian = self.condensed_matrix.calculate_rw_laplacian()
#        eigvectors = 

        clustering = SpectralClustering(k, mode ='amg')#,affinity='precomputed')
        clustering.fit(self.affinity_matrix)          
        algorithm_details = "Spectral Clustering (k = " +str(k)+")"
        clusters = gen_clusters_from_class_list(clustering.labels_)
        return Clustering(clusters,details = algorithm_details)
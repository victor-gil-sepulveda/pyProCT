'''
Created on 20/02/2012

@author: victor
'''

import validation.datasets as data
import scipy.spatial.distance as distance
from pyRMSD.condensedMatrix import CondensedMatrix
from pyproclust.algorithms.gromos.gromosAlgorithm import GromosAlgorithm
from pyproclust.algorithms.dbscan.dbscanAlgorithm import DBSCANAlgorithm
from pyproclust.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproclust.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
from pyproclust.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm
from validation.validationTools import params_to_string, dataset_loading_2D,\
    show_2D_dataset_clusters

"""
Script for visual validation of algorithms.
"""
def build_algorithms(matrix):
    algorithms = {
                  "GROMOS":GromosAlgorithm(matrix),
                  "DBSCAN":DBSCANAlgorithm(matrix),
                  "K-Medoids":KMedoidsAlgorithm(matrix),
                  "Hierarchical":HierarchicalClusteringAlgorithm(matrix),
                  "Spectral":SpectralClusteringAlgorithm(matrix)
                  }
    
    return algorithms

def generate_params_for_alg_and_dataset():
    
    params = {
              "GROMOS":{},
              "DBSCAN":{},
              "K-Medoids":{},
              "Hierarchical":{},
              "Spectral":{}
              }
    for algorithm_name in params.keys():
        for dataset_name in data.number_of_clusters:
            params[algorithm_name][dataset_name] = []
    
    # Params for spectral and K-medoids
    for dataset_name in data.number_of_clusters:
        num_clusters_list = data.number_of_clusters[dataset_name]
        print num_clusters_list
        for num_clusters in num_clusters_list:
            params["K-Medoids"][dataset_name].append({"k":num_clusters, "seeding_type":"EQUIDISTANT"})
            params["Spectral"][dataset_name].append({"k":num_clusters})
            
    # Params for GROMOS
    for dataset_name in data.number_of_clusters:
        params["GROMOS"][dataset_name].append({'cutoff':7.0})
        params["GROMOS"][dataset_name].append({'cutoff':5.0})
        params["GROMOS"][dataset_name].append({'cutoff':3.0})
        params["Hierarchical"][dataset_name].append({'cutoff':1.1523})
        
    return params

if __name__ == '__main__':
    params_for_alg_and_dataset = generate_params_for_alg_and_dataset()
    for dataset_name in data.all_datasets:
        
        dataset = data.all_datasets[dataset_name]
        
        # Creating the matrix
        observations = dataset_loading_2D(dataset)
        total_number_of_observations = len(observations)
        condensed_matrix_data = distance.pdist(observations)
        condensed_matrix = CondensedMatrix(condensed_matrix_data)
        
        algorithms = build_algorithms(condensed_matrix)
        for algorithm_name in algorithms:
            algorithm = algorithms[algorithm_name]
            params_list = params_for_alg_and_dataset[algorithm_name][dataset_name]
            for params in params_list:
                clustering = algorithm.perform_clustering(params)
                image_name = algorithm_name+"_"+dataset_name+"_"+params_to_string(params)
                show_2D_dataset_clusters(observations, 20, clustering, 20).save("clustering_images/%s.jpg"%image_name,"JPEG")

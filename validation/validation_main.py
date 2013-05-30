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
    show_2D_dataset_clusters, generate_similarity_network
from validation.datasets import sigma_sq
import numpy
from pyproclust.tools.scriptTools import create_directory

"""
Script for visual validation of algorithms.
"""
def get_algorithms():
    algorithms = {
                  "GROMOS":GromosAlgorithm,
                  "DBSCAN":DBSCANAlgorithm,
                  "K-Medoids":KMedoidsAlgorithm,
                  "Hierarchical":HierarchicalClusteringAlgorithm,
                  "Spectral":SpectralClusteringAlgorithm
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
        for num_clusters in num_clusters_list:
            params["K-Medoids"][dataset_name].append({"k":num_clusters, "seeding_type":"EQUIDISTANT"})
            params["Spectral"][dataset_name].append({"k":num_clusters, "use_k_medoids": True})
            
    # Params for GROMOS
    for dataset_name in data.number_of_clusters:
        params["GROMOS"][dataset_name].append({'cutoff':7.0})
        params["GROMOS"][dataset_name].append({'cutoff':5.0})
        params["GROMOS"][dataset_name].append({'cutoff':3.0})
        params["Hierarchical"][dataset_name].append({'cutoff':1.1523})
    
    params["DBSCAN"] = data.DBSCAN_params_sq
    
    return params

if __name__ == '__main__':
    
    params_for_alg_and_dataset = generate_params_for_alg_and_dataset()
    condensed_matrices = {}
    all_observations = {}
    create_directory("./clustering_images")
    for dataset_name in data.all_datasets:
        dataset = data.all_datasets[dataset_name]
        # Creating the matrix
        observations = dataset_loading_2D(dataset)
        all_observations[dataset_name] = observations
        condensed_matrix_data = distance.pdist(observations)
        condensed_matrix = CondensedMatrix(condensed_matrix_data)
        condensed_matrices[dataset_name] = condensed_matrix
        print "Matrix for %s:"%dataset_name
        print "-----------------------"
        print "Max dist. = ",condensed_matrix.calculateMax() 
        print "Min dist. = ",condensed_matrix.calculateMin()
        print "Mean dist. = ",condensed_matrix.calculateMean()
        print "Variance = ",condensed_matrix.calculateVariance()
        print "-----------------------\n"
    
    for dataset_name in data.all_datasets:
        print dataset_name
        observations = all_observations[dataset_name] 
        condensed_matrix = condensed_matrices[dataset_name]
        algorithms = get_algorithms()
        ###
        # And use Spectral.W to draw the adjacency graph in order to evaluate its correctness
        ###
        generate_similarity_network(algorithms["Spectral"].W,
                                    observations,30,20,False).save("clustering_images/%s_spectral_network.jpg"%dataset_name,
                                         "JPEG")
        ###
        for algorithm_name in algorithms:
            if(algorithm_name == "Spectral"):
                algorithm = algorithms[algorithm_name](condensed_matrix)
            else:
                algorithm = algorithms[algorithm_name](condensed_matrix)
            params_list = params_for_alg_and_dataset[algorithm_name][dataset_name]
            for params in params_list:
                clustering = algorithm.perform_clustering(params)
                image_name = algorithm_name+"_"+dataset_name+"_"+params_to_string(params)
                show_2D_dataset_clusters(observations, 20, clustering, 20).save("clustering_images/%s.jpg"%image_name,"JPEG")
    
    print
    print "Done"
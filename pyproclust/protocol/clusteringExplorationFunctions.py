'''
Created on 19/09/2012

@author: victor
'''
from pyproclust.algorithms.dbscan.dbscanTools import dbscan_param_space_search
from pyproclust.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproclust.algorithms.random.RandomAlgorithm import RandomClusteringAlgorithm
from pyproclust.algorithms.gromos.gromosAlgorithm import GromosAlgorithm
from pyproclust.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
from pyproclust.algorithms.hierarchical.hierarchicalTools import get_clusters_with_dicotomic_search
import pyproclust.tools.commonTools as common
from pyproclust.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm
from pyproclust.algorithms.dbscan.dbscanAlgorithm import DBSCANAlgorithm

def do_clustering_exploration(protocol_params, processManager,condensed_distance_matrix, max_dist, mean_dist,clustering_directory, htmlReport):
    
    if protocol_params.use_spectral:
        spawn_spectral(protocol_params, condensed_distance_matrix, processManager,clustering_directory)
    if protocol_params.use_gromos:
        spawn_gromos(protocol_params, condensed_distance_matrix,processManager,clustering_directory)
    if protocol_params.use_random:
        spawn_random(protocol_params, condensed_distance_matrix,processManager,clustering_directory)
    if protocol_params.use_hierarchical:
        spawn_hierarchical(protocol_params, condensed_distance_matrix, max_dist, processManager,clustering_directory) 
    if protocol_params.use_dbscan:
        spawn_dbscan(protocol_params, condensed_distance_matrix,processManager,clustering_directory)
    if protocol_params.use_kmedoids:
        spawn_kmedoids(protocol_params, condensed_distance_matrix, mean_dist, processManager,clustering_directory)
    processManager.consume() # Will wait until all the processes finish

def spawn_hierarchical(protocol_params, condensed_distance_matrix, max_dist, process_manager,clustering_directory):
    ###################
    # Hierarchical
    ###################
    common.print_and_flush("Spawning hierarchical clustering.\n")
    # Calculate cutoffs for hierarchical if necessary
    if len(protocol_params.hierarchical_cutoff_list) == 0: 
        spawn_hierarchical_cutoff_calculation(protocol_params,condensed_distance_matrix, max_dist, process_manager,clustering_directory)
    else:
        spawn_hierarchical_with_cutoffs(protocol_params, condensed_distance_matrix, max_dist, process_manager,clustering_directory)

def spawn_hierarchical_cutoff_calculation(protocol_params,condensed_distance_matrix, max_dist, process_manager,clustering_directory):
    process_name = "Hierarchical_Cutoff_Calculation"
    description = "Hierarchical Clustering Cutoff Calculation"
    function_kwargs={"protocol_params":protocol_params,"condensed_distance_matrix":condensed_distance_matrix,"max_dist":max_dist,"clustering_directory":clustering_directory}
    process_manager.add_process_internally(process_name,description,spawn_hierarchical_cutoff_calculation_function,function_kwargs,dependency = [])

def spawn_hierarchical_cutoff_calculation_function(protocol_params,condensed_distance_matrix,max_dist,clustering_directory):
    common.print_and_flush("Spawning cutoff calculation.\n")
    hierarchicalAlgorithm = HierarchicalClusteringAlgorithm(condensed_distance_matrix)
    common.print_and_flush("Calculating cutoffs for hierarchical...")
    clusters = get_clusters_with_dicotomic_search(condensed_distance_matrix,hierarchicalAlgorithm,0.,max_dist,protocol_params.min_clusters,protocol_params.max_clusters,protocol_params.hierarchical_cutoff_refinement_value)
    common.print_and_flush(" Done\n")
    
    common.print_and_flush("Storing clusters...")
    cid = 0
    for k in clusters.keys():
        clustering = clusters[k][1]
        complete_path_with_name = clustering_directory+"/"+"hie_"+str(cid)+".bin"
        clustering.save_to_disk(complete_path_with_name)
        cid = cid + 1
    common.print_and_flush(" Done\n")

def spawn_hierarchical_with_cutoffs(protocol_params, condensed_distance_matrix,process_manager,clustering_directory):
    hierarchicalAlgorithm = HierarchicalClusteringAlgorithm(condensed_distance_matrix)
    for c in protocol_params.hierarchical_cutoff_list:
        process_name = "Hierarchical_"+str(c)
        description = "Hierarchical clustering calculation with cutoff = "+str(c)
        function_kwargs={"algorithm":hierarchicalAlgorithm,"clustering_id":process_manager.next_process_id(),"algorithm_kwargs":{"cutoff":c},"directory":clustering_directory}
        process_manager.add_process_internally(process_name,description,run_algorithm,function_kwargs,dependency = [])

def spawn_gromos(protocol_params, condensed_distance_matrix,process_manager,clustering_directory):
    ###################
    # Gromos
    ###################
    common.print_and_flush("Spawning gromos clustering.\n")
    alg = GromosAlgorithm(condensed_distance_matrix)
    for c in protocol_params.gromos_cutoff_list:
        process_name = "Gromos_"+str(c)
        description = "GROMOS Algorithm Calculation. Cutoff:"+str(c)
        function_kwargs={"algorithm":alg,"clustering_id":process_manager.next_process_id(),"algorithm_kwargs":{"cutoff":c},"directory":clustering_directory}
        process_manager.add_process_internally(process_name,description,run_algorithm,function_kwargs,dependency = [])

def spawn_random( protocol_params, condensed_distance_matrix,process_manager,clustering_directory):
    ###################
    # Random
    ###################
    common.print_and_flush("Spawning random clustering.\n")
    alg = RandomClusteringAlgorithm(condensed_distance_matrix)
    for i in range(protocol_params.max_random_clusters):
        process_name = "Random_" + str(i)
        description = "Random Algorithm Calculation "+str(i)
        function_kwargs={"algorithm":alg, "clustering_id":process_manager.next_process_id(), "algorithm_kwargs":{"max_num_of_clusters":protocol_params.max_clusters},"directory":clustering_directory}
        process_manager.add_process_internally(process_name,description,run_algorithm,function_kwargs,dependency = [])

def spawn_kmedoids(protocol_params, condensed_distance_matrix,starting_gromos_cutoff,process_manager,clustering_directory):
    ###################
    # KMEDOIDS
    ###################
    common.print_and_flush("Spawning k-medoids clustering.\n")
    base_cutoff = starting_gromos_cutoff
    alg = KMedoidsAlgorithm(condensed_distance_matrix)
    for k in range(protocol_params.min_clusters, protocol_params.max_clusters+1, protocol_params.kmedoids_step):
        process_name = "KMedoids_" + str(k)
        description = "K-Medoids calculation for k = "+str( k)+" and base seeding cutoff "+str(base_cutoff)
        function_kwargs={"algorithm":alg, "clustering_id":process_manager.next_process_id(), "algorithm_kwargs":{"k":k,"seeding_max_cutoff":base_cutoff},"directory":clustering_directory}
        process_manager.add_process_internally(process_name,description,run_algorithm,function_kwargs,dependency = [])
        base_cutoff -= 0.1

def spawn_dbscan(protocol_params, condensed_distance_matrix,process_manager,clustering_directory):
    ###################
    # DBSCAN
    ###################
    dbscan_param_pairs = []
    if len(protocol_params.dbscan_param_pairs)==0:
        common.print_and_flush( "Choosing better parameters for DBSCAN ... ")
        dbscan_param_pairs = dbscan_param_space_search(condensed_distance_matrix.row_length, protocol_params.max_noise, condensed_distance_matrix)
        common.print_and_flush( " Done\n")
    else:
        dbscan_param_pairs = protocol_params.dbscan_param_pairs
    
    common.print_and_flush("Spawning DBSCAN clustering.\n")
    alg = DBSCANAlgorithm(condensed_distance_matrix)
    for pair in dbscan_param_pairs:
        process_name = "DBSCAN_minpts_" + str(pair[0])+"_cutoff_"+str(pair[1])
        description = "DBSCAN calculation with minpts = " + str(pair[0])+" and  cutoff = "+str(pair[1])
        function_kwargs={"algorithm":alg, "clustering_id":process_manager.next_process_id(), "algorithm_kwargs":{"eps": pair[1], "minpts":pair[0]},"directory":clustering_directory}
        process_manager.add_process_internally(process_name,description,run_algorithm,function_kwargs,dependency = [])

def spawn_spectral(protocol_params, condensed_distance_matrix, processManager,clustering_directory):
    ###################
    # SPECTRAL CLUSTERING
    ###################
    algorithm = SpectralClusteringAlgorithm(condensed_distance_matrix)
    for k in range(protocol_params.min_clusters,protocol_params.max_clusters+1,protocol_params.spectral_clustering_step):
        process_name = "Spectral_clustering_k_" + str(k)
        description = process_name
        function_kwargs={"algorithm":algorithm, "clustering_id":processManager.next_process_id(), "algorithm_kwargs":{"k": k},"directory":clustering_directory}
        if k != protocol_params.min_clusters: # If it's not the first
            dependency = ["Spectral_clustering_k_" + str(k-protocol_params.spectral_clustering_step)]
        else:
            dependency = []
        processManager.add_process_internally(process_name,description,run_algorithm,function_kwargs,dependency)

def run_algorithm(algorithm, algorithm_kwargs,clustering_id,directory):
    clustering = algorithm.perform_clustering(algorithm_kwargs)
    clustering.save_to_disk(directory+"/"+str(clustering_id)+".bin")

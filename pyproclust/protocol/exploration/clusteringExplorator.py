'''
Created on 05/02/2013

@author: victor
'''
from pyproclust.gui.observer.Observable import Observable
from pyproclust.protocol.processPool import ProcessPool
from pyproclust.protocol.serialProcessPool import SerialProcessPool
from pyproclust.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm
from pyproclust.algorithms.dbscan.dbscanAlgorithm import DBSCANAlgorithm
from pyproclust.algorithms.gromos.gromosAlgorithm import GromosAlgorithm
from pyproclust.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproclust.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm
from pyproclust.algorithms.random.RandomAlgorithm import RandomClusteringAlgorithm
from pyproclust.clustering.clustering import Clustering
import os


def run_algorithm(algorithm, algorithm_kwargs, clustering_id, directory):
    clustering = algorithm.perform_clustering(algorithm_kwargs)
    clustering.save_to_disk(directory+"/"+str(clustering_id)+".bin")

class ClusteringExplorator(Observable):

    @classmethod
    def get_available_algorithms(cls):
        return { 
                "spectral": SpectralClusteringAlgorithm,
                "dbscan":  DBSCANAlgorithm,
                "gromos": GromosAlgorithm,
                "kmedoids": KMedoidsAlgorithm,
                "random": RandomClusteringAlgorithm,
                "hierarchical": HierarchicalClusteringAlgorithm
        }
    
    def __init__(self, parameters, matrix_handler, workspace_handler, scheduler_type, observer):
        super(ClusteringExplorator,self).__init__(observer)
        
        self.matrix_handler = matrix_handler
        self.workspace_handler = workspace_handler
        self.parameters = parameters
        
        if scheduler_type == "Process/Parallel":
            self.scheduler = ProcessPool(parameters["control"]["number_of_processors"], 
                                         parameters["control"]["algorithm_scheduler_sleep_time"])
        elif scheduler_type == "Serial":
            self.scheduler = SerialProcessPool()
        else:
            print "[ERROR][ClusteringExplorator::__init__] Not supported scheduler_type ( %s )"%(scheduler_type)
            exit()
        
        self.current_clustering_id = 0
    
    def get_used_algorithms(self, parameters):
        names = []
        for algorithm_name in parameters["algorithms"]:
            if parameters["algorithms"][algorithm_name]["use"]:
                names.append(algorithm_name)
        return names
    
    def run(self):
        used_algorithms = self.get_used_algorithms(self.parameters)
        
        # Generate all clustering + info structures
        clusterings_info = {}
        for algorithm_type in used_algorithms:
            clusterings_info = dict(clusterings_info.items() + self.schedule_algorithms( algorithm_type).items())
        
        # Wait until all processes have finished
        self.scheduler.consume()
        
        # Load clusterings and put them inside the structure
        clustering_plus_files = Clustering.load_all_from_directory(self.workspace_handler["clusterings"])
        for clustering, filename_with_extension in clustering_plus_files:
            clustering_id = os.path.split(filename_with_extension)[1].split(".")[0]
            clusterings_info[clustering_id]["clustering"] = clustering
        
        return clusterings_info
    
    def schedule_algorithms(self, algorithm_type):
        
        algorithm_data = self.parameters["algorithms"][algorithm_type]
        
        # The algorithm we are going to use
        algorithm = self.build_algorithm(algorithm_type)
        
        # A list with all the parameters for diverse runs
        algorithm_run_params = algorithm_data["parameters"]
        
        # If not parameters were given we have to get the better ones
        clusterings_info = {}
        clusterings_have_been_recalculated  = False
        
        if algorithm_data["auto"]:
            clusterings_info, clusterings_have_been_recalculated =  self.generate_parameters_for_algorithm(algorithm_type)
        else:
            clusterings_info =  self.generate_clustering_info(algorithm_type, algorithm_run_params)
        
        # Sometimes getting the best parameters imply getting the clusterings itselves
        if not clusterings_have_been_recalculated:
            for clustering_id in clusterings_info:
                one_clustering_info = clusterings_info[clustering_id]
                self.scheduler.add_process_internally(clustering_id,
                                                      "Generation of clustering with %s algorithm with id %s"%(
                                                                                                one_clustering_info["type"],
                                                                                                clustering_id
                                                                                                ),
                                                      run_algorithm,
                                                                  {
                                                                   "algorithm":algorithm, 
                                                                   "clustering_id":clustering_id, 
                                                                   "algorithm_kwargs":one_clustering_info["parameters"],
                                                                   "directory":self.workspace_handler["clusterings"]
                                                                   },
                                                      dependencies = [])
        return clusterings_info
    
    def generate_clustering_info(self, algorithm_type, clustering_parameters):
        clustering_info = {}
        for running_parameters in clustering_parameters:
            clustering_id = "clustering_%04d"%(self.current_clustering_id)
            self.current_clustering_id += 1
            clustering_info[clustering_id] = {
                                                "type":algorithm_type,
                                                "clustering": None,
                                                "parameters": running_parameters
            }
        return clustering_info
     
    def generate_parameters_for_algorithm(self, algorithm_type):
        if not algorithm_type in ["spectral","dbscan","gromos","kmedoids","random","hierarchical"] :
            print "[ERROR][ClusteringExplorator::build_algorithms] Not known algorithm type( %s )"%(algorithm_type)
            exit()
        
    
    def build_algorithm(self, algorithm_type):
        distance_matrix = self.matrix_handler.distance_matrix
        
        # We need to set number of clusters for performance and to get sigma
        algorithm_execution_parameters = {"max_clusters":self.parameters["evaluation"]["maximum_clusters"]}
        try:
            algorithm_execution_parameters["sigma_sq"] = self.parameters["algorithms"]["spectral"]["sigma"]
        except KeyError:
            pass
        
        if algorithm_type in ["spectral","dbscan","gromos","kmedoids","random","hierarchical"] :
            return ClusteringExplorator.get_available_algorithms()[algorithm_type](distance_matrix, **algorithm_execution_parameters)
        
        else:
            print "[ERROR][ClusteringExplorator::build_algorithms] Not known algorithm type( %s )"%(algorithm_type)
            exit()
    
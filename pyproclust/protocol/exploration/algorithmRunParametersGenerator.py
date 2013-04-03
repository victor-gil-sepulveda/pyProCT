'''
Created on 06/02/2013

@author: victor
'''
import pyproclust.algorithms.dbscan.dbscanTools as dbscanTools
import pyproclust.algorithms.hierarchical.hierarchicalTools as  hierarchicalTools
from pyproclust.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm

class AlgorithmRunParametersGenerator(object):

    CUTOFF_STEP = 0.25
    CLUSTERING_SIZE_STEP = 2
    HIERARCHICAL_REFINEMENT_VALUE = 200
    
    def __init__(self, parameters, matrix_handler):
        """
        Class creator.
        
        @param parameters: Script parameters.
        
        @param distance_matrix: The distance matrix we are using.
        """
        self.distance_matrix = matrix_handler.distance_matrix
        self.parameters = parameters
    
    @classmethod
    def get_base_parameters(cls):
        """
        Defines the base parameters needed for each of the algorithms. Each parameter created will be based
        on one of those and must not have more keys than these.
        
        @return: A dictionary with the base parameters needed for each of the algorithms. 
        """
        return { 
                 "gromos":{
                        "cutoff": None
                  },
                  "spectral":{
                         "k": None,
                         "use_k_medoids": None
                  },
                  "dbscan":{
                         "eps": None,
                         "minpts": None
                  },
                  "hierarchical":{
                         "cutoff": None,
                         "method": None
                  },
                  "kmedoids":{
                         "k": None,
                         "seeding_type": None,
                         "seeding_max_cutoff": None
                  },
                  "random":{
                        "num_clusters": None,
                        
                  }
        }
    
    def get_parameters_for_type(self, algorithm_type):
        """
        Tries to get some parameters to do a good exploration using one algorithm.
        
        @param algorithm_type: A valid algorithm type for which to get the parameters.
        
        @return: A tuple consisting on the parameters we are going to use and, if the clusterings were
        generated to get the clusterings, then those clusterings are also returned for the sake of performance.
        In this case the algorithm parameters in position 'i' correspond to the created cluster in position 'i'. 
        """
        
        if algorithm_type == "gromos":
            return self.get_gromos_parameters()
                
        elif algorithm_type == "spectral":
            return self.get_spectral_parameters()
                
        elif algorithm_type == "kmedoids":
            return self.get_kmedoids_parameters()
        
        elif algorithm_type == "hierarchical":
            return self.get_hierarchical_parameters()
                
        elif algorithm_type == "dbscan":
            return self.get_dbscan_parameters()
                
        elif algorithm_type == "random":
            return self.get_random_parameters()
       
        return None

    def get_gromos_parameters(self):
        """
        This function creates some parameters to be used with Gromos. 
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []
        cutoff = AlgorithmRunParametersGenerator.CUTOFF_STEP
        while cutoff < self.distance_matrix.calculateMean():
            run_parameter = AlgorithmRunParametersGenerator.get_base_parameters()["gromos"]
            run_parameter["cutoff"] = cutoff
            cutoff += AlgorithmRunParametersGenerator.CUTOFF_STEP
            run_parameters.append(run_parameter)
        
        return run_parameters, []
        
    def get_spectral_parameters(self):
        """
        This function creates some parameters to be used with the spectral clustering algorithm. 
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []
        max_clusters = self.parameters["evaluation"]["maximum_clusters"]
        min_clusters = self.parameters["evaluation"]["minimum_clusters"]
        sizes = range(min_clusters,max_clusters,AlgorithmRunParametersGenerator.CLUSTERING_SIZE_STEP)
        for one_size in sizes:
            run_parameter = AlgorithmRunParametersGenerator.get_base_parameters()["spectral"]
            run_parameter["k"]  = one_size
            run_parameter["use_k_medoids"] = True
            run_parameters.append(run_parameter)
    
        return run_parameters, []
    
    def get_kmedoids_parameters(self):
        """
        This function creates some parameters to be used with K-Medoids. 
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []
        max_clusters = self.parameters["evaluation"]["maximum_clusters"]
        min_clusters = self.parameters["evaluation"]["minimum_clusters"]
        sizes = range(min_clusters,max_clusters,AlgorithmRunParametersGenerator.CLUSTERING_SIZE_STEP)
        for one_size in sizes:
            run_parameter = AlgorithmRunParametersGenerator.get_base_parameters()["kmedoids"]
            run_parameter["k"]  = one_size
            run_parameter["seeding_type"] = "GROMOS"
            run_parameter["seeding_max_cutoff"] = self.distance_matrix.calculateMean()
            run_parameters.append(run_parameter)
        
        return run_parameters, []
    
    def get_hierarchical_parameters(self):
        """
        This function creates some parameters to be used with the Hierarchical algorithm. 
        @return: A tuple with the generated parameters and a list with its corresponding clusterings 
        (in this case parameters and clusterings are obtained at the same time).
        """
        run_parameters = []
        max_clusters = self.parameters["evaluation"]["maximum_clusters"]
        min_clusters = self.parameters["evaluation"]["minimum_clusters"]
        hierarchicalAlgorithm = HierarchicalClusteringAlgorithm(self.distance_matrix)
        clusters_and_cutoff = hierarchicalTools.get_clusters_with_dicotomic_search(self.distance_matrix,
                                                                        hierarchicalAlgorithm,
                                                                        0.,
                                                                        self.distance_matrix.calculateMax(),
                                                                        min_clusters,
                                                                        max_clusters,
                                                                        AlgorithmRunParametersGenerator.HIERARCHICAL_REFINEMENT_VALUE)
        clusterings = []
        cutoffs = []
        for numclusters in clusters_and_cutoff:
            clusterings.append(clusters_and_cutoff[numclusters][1])
            cutoffs.append(clusters_and_cutoff[numclusters][0])
#         (clusterings, cutoffs) = zip(*clusters_and_cutoff)
        
        for cutoff in cutoffs:
            run_parameter = AlgorithmRunParametersGenerator.get_base_parameters()["hierarchical"]
            run_parameter["method"] = 'complete'
            run_parameter["cutoff"] = cutoff
            run_parameters.append(run_parameter)
        
        return run_parameters, clusterings
    
    def get_dbscan_parameters(self):
        """
        This function creates some parameters to be used with DBScan. 
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []
        # (minpts, eps tuples)
        dbscan_param_pairs = dbscanTools.dbscan_param_space_search(self.parameters["evaluation"]["maximum_noise"], 
                                                                   self.distance_matrix)
        for (minpts, eps) in dbscan_param_pairs:
            run_parameter = AlgorithmRunParametersGenerator.get_base_parameters()["dbscan"]
            run_parameter["minpts"] = minpts
            run_parameter["eps"] = eps
            run_parameters.append(run_parameter)
        
        return run_parameters, []
    
    def get_random_parameters(self):
        """
        This function creates some parameters to be used with the random algorithm. 
        @return: A tuple with the generated parameters and an empty list corresponding to the clusterings.
        """
        run_parameters = []
        max_clusters = self.parameters["evaluation"]["maximum_clusters"]
        min_clusters = self.parameters["evaluation"]["minimum_clusters"]
        sizes = range(min_clusters,max_clusters,AlgorithmRunParametersGenerator.CLUSTERING_SIZE_STEP)
        for one_size in sizes:
            run_parameter = AlgorithmRunParametersGenerator.get_base_parameters()["random"]
            run_parameter["num_clusters"]  = one_size
            run_parameters.append(run_parameter)
        
        return run_parameters, []
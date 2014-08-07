"""
Created on 05/02/2013

@author: victor
"""
from pyproct.driver.observer.observable import Observable
from pyproct.clustering.algorithms.spectral.spectralClusteringAlgorithm import SpectralClusteringAlgorithm
from pyproct.clustering.algorithms.dbscan.dbscanAlgorithm import DBSCANAlgorithm
from pyproct.clustering.algorithms.gromos.gromosAlgorithm import GromosAlgorithm
from pyproct.clustering.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproct.clustering.algorithms.random.RandomAlgorithm import RandomClusteringAlgorithm
from pyproct.clustering.algorithms.hierarchical.hierarchicalAlgorithm import HierarchicalClusteringAlgorithm

def run_algorithm(algorithm, algorithm_kwargs, clustering_id):
    """
    This function launches an execution of one clustering algorithm with its parameters. Used mainly to be
    scheduled.

    @param algorithm: Instance of a clustering algorithm.

    @param algorithm_kwargs: The parameters needed by the algorithm above to run.

    @param clustering_id: An id used to define the resulting clustering.
    """
    clustering = algorithm.perform_clustering(algorithm_kwargs)
    return (clustering_id, clustering)

class ClusteringExplorer(Observable):

    @classmethod
    def get_clustering_algorithm_class(cls):
        """
        This function returns a dictionary that links a clustering algorithm type with its class creator.
        @return: The aforementioned dictionary.
        """
        return {
                "spectral": SpectralClusteringAlgorithm,
                "dbscan":  DBSCANAlgorithm,
                "gromos": GromosAlgorithm,
                "kmedoids": KMedoidsAlgorithm,
                "random": RandomClusteringAlgorithm,
                "hierarchical": HierarchicalClusteringAlgorithm
        }

    def __init__(self, parameters, matrix_handler, workspace_handler, scheduler, parameters_generator, observer = None):
        """
        Class creator.

        @param parameters: Script parameters.

        @param matrix_handler: The matrix handler MatrixHandler instance) containing the distance matrix.

        @param workspace_handler: The workspace handler containing paths.

        @param scheduler: An instance of a Scheduler like object

        @param parameters_generator: An instance of AlgorithmRunParametersGenerator, in charge of generating automatically the
        parameters needed for the exploration if those are not given.

        @param observer: The observer object for this Observable.
        """
        super(ClusteringExplorer,self).__init__(observer)

        self.matrix_handler = matrix_handler
        self.workspace_handler = workspace_handler
        self.clustering_parameters = parameters["clustering"]
        self.evaluation_parameters = parameters["clustering"]["evaluation"]
        self.current_clustering_id = 0
        self.parameters_generator = parameters_generator
        self.scheduler = scheduler

    def run(self):
        """
        Executes the whole exploration pipeline:
            - Generates the parameters structures
            - Executes the algorithms for all different parameters
            - Loads the results

        @return: A dictionary 'clustering_info' structures indexed by clustering ID. Each of these structures
        contains one generated clustering as well as the algorithm type and parameters used to get it.
        """
        used_algorithms = self.clustering_parameters["algorithms"].keys()

        # Generate all clustering + info structures
        clusterings_info = {}
        for algorithm_type in used_algorithms:
            clusterings_info = dict(clusterings_info.items() + self.schedule_algorithm(algorithm_type).items())# append elements to a dict

        # Wait until all processes have finished
        clusterings = self.scheduler.run()

        # Put clusterings inside the structure
        for clustering_id, clustering in clusterings:
            clusterings_info[clustering_id]["clustering"] = clustering

        return clusterings_info

    def schedule_algorithm(self, algorithm_type):
        """
        Structures all the info needed for an algorithm+parameters execution and pushes it into the scheduling queue.

        @param algorithm_type: The algorithm type of the clustering algorithm we are working with.

        @return: The 'clustering info' structure as defined by  AlgorithmRunParametersGenerator::get_parameters_for_type.
        """

        algorithm_data = self.clustering_parameters["algorithms"][algorithm_type]

        # The algorithm we are going to use
        algorithm = self.build_algorithm(algorithm_type)

        # If not parameters were given we have to get the better ones
        clusterings = []

        auto_parameter_generation = True if not "parameters" in algorithm_data else False

        if auto_parameter_generation:
            print "Generating params for", algorithm_type
            algorithm_run_params, clusterings =  self.parameters_generator.get_parameters_for_type(algorithm_type)
        else:
            # A list with all the parameters for diverse runs
            algorithm_run_params = algorithm_data["parameters"]

        clusterings_info =  self.generate_clustering_info(algorithm_type, algorithm_run_params, clusterings)

        # Sometimes getting the best parameters imply getting the clusterings themselves
        if clusterings == []:
            for clustering_id in clusterings_info:
                one_clustering_info = clusterings_info[clustering_id]
                self.scheduler.add_task( task_name = clustering_id,
                                         description = "Generation of clustering with %s algorithm and id %s"%(
                                                                                    one_clustering_info["type"],
                                                                                    clustering_id
                                                                                    ),
                                          target_function = run_algorithm,
                                          function_kwargs = {
                                                       "algorithm":algorithm,
                                                       "clustering_id":clustering_id,
                                                       "algorithm_kwargs":one_clustering_info["parameters"]
                                                       },
                                          dependencies = {})
        return clusterings_info

    def generate_clustering_info(self, algorithm_type, clustering_parameters, clusterings = []):
        """
        It builds the clustering_info structures by parsing the parameters.

        @param algorithm_type: The algorithm type of the clustering algorithm we are working with.

        @param clustering_parameters: The parameters we are going to try with this algorithm.

        @param clusterings: In the case that the parameter generation also created the clusterings, this argument
        will hold them. Clustering parameters and clusterings must be correlated so that 'clustering_parameters[i]' where
        the parameters used to get 'clustering[i]'.

        @return: A list of clustering_info structures.
        """
        clustering_info = {}
        for i, running_parameters in enumerate(clustering_parameters):

            clustering_id = "clustering_%04d"%(self.current_clustering_id)
            self.current_clustering_id += 1
            clustering_info[clustering_id] = {
                                                "type":algorithm_type,
                                                "clustering": None,
                                                "parameters": running_parameters
            }

            if clusterings != []:
                clustering_info[clustering_id]["clustering"] = clusterings[i]

        return clustering_info

    def build_algorithm(self, algorithm_type):
        """
        Creates an algorithm with type 'algorithm_type'.

        @param algorithm_type: The algorithm type.

        @return: An instance of an algorithm of type 'algorithm_type'.
        """
        distance_matrix = self.matrix_handler.distance_matrix
        algorithm_execution_parameters = {}
        if algorithm_type == "spectral":
            # We need to set number of clusters for performance and we get sigma if defined
            algorithm_execution_parameters["max_clusters"] = self.evaluation_parameters["maximum_clusters"]
            if "sigma" in self.clustering_parameters["algorithms"]["spectral"]:
                algorithm_execution_parameters["sigma_sq"] = self.clustering_parameters["algorithms"]["spectral"]["sigma"]
            # else it calculates its own sigma

        if algorithm_type in ["spectral","dbscan","gromos","kmedoids","random","hierarchical"] :
            return ClusteringExplorer.get_clustering_algorithm_class()[algorithm_type](distance_matrix, **algorithm_execution_parameters)
        else:
            print "[ERROR][ClusteringExplorer::build_algorithms] Not known algorithm type ( %s )"%(algorithm_type)
            self.notify("SHUTDOWN", "Not known algorithm type ( %s )"%(algorithm_type))
            exit()

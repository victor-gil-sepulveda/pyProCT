"""
Created on 11/02/2013

@author: victor
"""
from pyproct.tools.matrixTools import get_submatrix
from pyproct.clustering.cluster import Cluster
from pyproct.clustering.clustering import Clustering
import pyproct.driver.scheduling.tools as scheduling_tools
from pyproct.clustering.selection.bestClusteringSelector import BestClusteringSelector
from pyproct.driver.observer.observable import Observable
from pyproct.clustering.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproct.clustering.evaluation.analysis.analysisRunner import AnalysisRunner
from pyproct.clustering.evaluation.analysis.analysisPopulator import AnalysisPopulator

class Refiner(Observable):
    """
    Given a clustering, it tries to further separate its clusters. This way we expect to reduce the scaling problem.
    """
    KMedoidsAlgorithmClass = KMedoidsAlgorithm

    def __init__(self, matrixHandler,  trajectoryHandler, clustering_parameters, refinement_parameters, observer):
        """
        """
        super(Refiner,self).__init__(observer)
        self.matrixHandler = matrixHandler
        self.trajectoryHandler = trajectoryHandler
        self.clustering_parameters = clustering_parameters
        self.refinement_paramenters = refinement_parameters


    @classmethod
    def redefine_cluster_with_map(cls, initial_cluster, new_cluster):
        """
        It renames the elements of a cluster using a map array.For instance, if our dataset has 5
        elements and the cluster to remap contains element 1, 3 and 4, then the elements map will
        remap 0->1, 1->3 and 2->4. This mapping is implemented as an array so map[0] = 1, map[1] = 3
        and map[2] = 4.
        """
        elements_map = initial_cluster.all_elements # The map are the elements themselves
        remapped_new_cluster_elems = []
        for element in new_cluster.all_elements:
            remapped_new_cluster_elems.append(elements_map[element])

        return Cluster(None, remapped_new_cluster_elems)

    @classmethod
    def repartition_with_kmedoids(cls, initial_cluster, k, submatrix):
        partitioned_clustering = cls.KMedoidsAlgorithmClass(submatrix).perform_clustering({
                                          "k": k,
                                          "seeding_type": "RANDOM",
                                          "tries":10})
        remapped_clusters = []
        for partitioned_cluster in partitioned_clustering.clusters:
            remapped_clusters.append(cls.redefine_cluster_with_map(initial_cluster, partitioned_cluster))
        return Clustering(remapped_clusters)

    def run (self, clustering):
        """
        Refine a clustering recursively using a k-means over each cluster.
        New clusters obtained from a cluster must have no noise and
        """
        max_partitions = self.refinement_parameters["max_partitions"]
        try_step = int(max(1, float(max_partitions) / self.refinement_parameters["tries_per_cluster"]))
        matrix = self.matrixHandler.distance_matrix

        new_clusters = []
        for cluster in clustering.clusters:
            base_id = cluster.id
            # The initial clustering is added to the list of new clusters.
            # With this 'trick' the initial cluster also enters the competition for the best clustering price.
            clusterings = {base_id:{"type":"refined_base",
                                    "clustering": Clustering([cluster]),
                                    "parameters": {}}}

            submatrix = get_submatrix(matrix, cluster.all_elements)

            # Proceed with some K Medoids partitions
            # TODO: Generate parameters with parameter generator
            for k in range(2,max_partitions,try_step):
                clustering = self.repartition_with_kmedoids(cluster, k, submatrix)
                clusterings["%s_%d"%(base_id,k)] = {"type":"refined",
                                                     "clustering": clustering,
                                                     "parameters": {"k":k}}

            # Evaluate all clusterings and pick the best one
            AnalysisRunner(scheduling_tools.build_scheduler(
                                                       self.clustering_parameters["clustering"]["control"],
                                                       self.observer),
                                          clusterings,
                                          AnalysisPopulator(self.matrixHandler,
                                                            self.trajectoryHandler,
                                                            self.clustering_parameters)).evaluate()

            best_clustering_id, all_scores = BestClusteringSelector(self.clustering_parameters).choose_best(clusterings)  # @UnusedVariable
            new_clusters.extend(clusterings[best_clustering_id]["clustering"].clusters)

        # Convert all new clusters in the new clustering
        return {"type":"refined_clustering",
                "clustering": Clustering(new_clusters),
                "parameters": self.refinement_parameters}





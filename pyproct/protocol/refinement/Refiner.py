'''
Created on 11/02/2013

@author: victor
'''
from pyproct.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm
from pyproct.tools.matrixTools import get_submatrix
from pyproct.clustering.cluster import Cluster


class Refiner():
    """
    Given a clustering, it tries to further separate its clusters. This way we expect to reduce the scaling problem.
    """

    def __init__(self):
        """
        """
        pass

    @classmethod
    def get_cluster_submatrix_and_map(cls, cluster, matrix):
        """
        Generates a condensed matrix from the original matrix and a map to recover element ids.
        For instance, if our dataset has 5 elements and the cluster to remap contains element 1,
        3 and 4, then the elements map will remap 0->1, 1->3 and 2->4. This mapping is implemented
        as an array so map[0] = 1, map[1] = 3 and map[2] = 4.

        """
        # Get cluster submatrix (like with compression)
        elements_map = []
        for e in cluster.all_elements:
            elements_map.append(e) # Position [0] for element x1 etc...
        submatrix = get_submatrix(matrix, cluster.all_elements)
        return submatrix, elements_map

    @classmethod
    def redefine_cluster_with_map(initial_cluster, new_cluster):
        """
        It renames the elements of a cluster using a map array.
        """
        elements_map = initial_cluster.all_elements # The map are the elements themselves
        remapped_new_cluster_elems = []
        for element in new_cluster.all_elements:
            remapped_new_cluster_elems.append(elements_map[element])

        return Cluster(None, remapped_new_cluster_elems)

    def run(self, clustering):
        """
        Refine a clustering recursively using a k-means over each cluster.
        New clusters obtained from a cluster must have no noise and
        """
#         clusters_at_most =
#         try_step =
#         max_random_tries =
#         new_clusters = []
#         for cluster in clustering.clusters:
#             # The initial clustering is added to the list of new clusters.
#             # With this 'trick' the initial cluster also enters the competition for the best clustering price.
#             clusterings = [Clustering([cluster])]
#
#             submatrix, elements_map = self.get_cluster_submatrix_and_map(cluster)
#
#             # Proceed with some K Medoids partitions
#             kmedoids = KMedoidsAlgorithm(submatrix)
#             for k in range(2,clusters_at_most,try_step):
#                 for k in range(max_random_tries):
#                     clusters = kmedoids.perform_clustering({
#                                                       "k": k,
#                                                       "seeding_type": "RANDOM"
#                     }).clusters
#                     clusterings.append(Clustering(redefine_clusters_with_map(clusters, elements_map)))
#
#             # Remove equal clusterings and those with noise
#
#
#             # Evaluate all clusterings
#             AnalysisRunner(scheduling_tools.build_scheduler(
#                                                        clustering_parameters["clustering"]["control"],
#                                                        self.observer),
#                                           clusterings,
#                                           AnalysisPopulator(matrixHandler,
#                                                             trajectoryHandler,
#                                                             clustering_parameters)).evaluate()
#
#             best_clustering_id, all_scores = BestClusteringSelector(clustering_parameters).choose_best(clusterings)

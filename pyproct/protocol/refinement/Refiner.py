'''
Created on 11/02/2013

@author: victor
'''
from pyproct.algorithms.kmedoids.kMedoidsAlgorithm import KMedoidsAlgorithm


class Refiner():

    def __init__(self):
        pass

    def run(self, clustering):
        """
        Refine a clustering recursively using a k-means over each cluster.
        New clusters obtained from a cluster must have no noise and
        """
        clusters_at_most =
        try_step =
        max_random_tries =
        new_clusters = []
        for cluster in clustering.clusters:
            # The initial clustering is added to the list of new clusters.
            # With this 'trick' we allow the best clustering to be the original one.
            clusterings = [Clustering([cluster])]

            # Get cluster submatrix (like with compression)
            elements_map = []
            for e in cluster.all_elements:
                elements_map.append(e) # Position [0] for element x1 etc...
            submatrix = recreate_matrix(matrix, cluster_size, elements_map)

            # Proceed with some K Medoids partitions
            kmedoids = KMedoidsAlgorithm(submatrix)
            for k in range(2,clusters_at_most,try_step):
                for k in range(max_random_tries):
                    clusters = kmedoids.perform_clustering({
                                                      "k": k,
                                                      "seeding_type": "RANDOM"
                    }).clusters
                    clusterings.append(Clustering(redefine_clusters_with_map(clusters, elements_map)))

            # Remove equal clusterings and those with noise


            # Evaluate all clusterings
            AnalysisRunner(scheduling_tools.build_scheduler(
                                                       clustering_parameters["clustering"]["control"],
                                                       self.observer),
                                          clusterings,
                                          AnalysisPopulator(matrixHandler,
                                                            trajectoryHandler,
                                                            clustering_parameters)).evaluate()

            best_clustering_id, all_scores = BestClusteringSelector(clustering_parameters).choose_best(clusterings)

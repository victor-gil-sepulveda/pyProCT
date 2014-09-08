"""
Created on 21/05/2012

@author: victor
"""
import pyproct.driver.scheduling.tools as scheduling_tools
from pyproct.clustering.selection.bestClusteringSelector import BestClusteringSelector
from pyproct.driver.observer.observable import Observable
from pyproct.clustering.protocol.exploration.clusteringExplorer import ClusteringExplorer
from pyproct.clustering.protocol.exploration.algorithmRunParametersGenerator import AlgorithmRunParametersGenerator
from pyproct.clustering.filtering.clusteringFilter import ClusteringFilter
from pyproct.clustering.evaluation.analysis.analysisRunner import AnalysisRunner
from pyproct.clustering.evaluation.analysis.analysisPopulator import AnalysisPopulator

class ClusteringProtocol(Observable):

    def __init__(self, timer, observer):
        super(ClusteringProtocol, self).__init__(observer)
        self.timer = timer

    def run(self, clustering_parameters, matrix_handler, data_handler, workspaceHandler):

        ############################
        # Clustering exploration
        ############################
        self.notify("Exploration Started", [])
        self.timer.start("Clustering Exploration")
        clusterings  = ClusteringExplorer(  clustering_parameters,
                                            matrix_handler,
                                            workspaceHandler,
                                            scheduling_tools.build_scheduler(clustering_parameters["global"]["control"],
                                                                             self.observer),
                                            AlgorithmRunParametersGenerator(clustering_parameters,
                                                                            matrix_handler),
                                            self.observer).run()

        self.notify("Clusterings Created", {"number_of_clusters":len(clusterings)})
        self.timer.stop("Clustering Exploration")

        ######################
        # First filtering
        ######################
        self.timer.start("Clustering Filtering")
        selected_clusterings, not_selected_clusterings = ClusteringFilter(clustering_parameters["clustering"]["evaluation"],
                                                                          matrix_handler).filter(clusterings)

        self.notify("Filter", {"selected":len(selected_clusterings.keys()),"not_selected":len(not_selected_clusterings.keys())})
        self.timer.stop("Clustering Filtering")

        if selected_clusterings == {}:
            return None

        ######################
        # Clustering scoring
        ######################
        self.timer.start("Evaluation")
        analyzer = AnalysisRunner(scheduling_tools.build_scheduler(
                                                       clustering_parameters["global"]["control"],
                                                       self.observer),
                                          selected_clusterings,
                                          AnalysisPopulator(matrix_handler,
                                                            data_handler,
                                                            clustering_parameters))

        analyzer.evaluate()
        self.timer.stop("Evaluation")

        ######################
        # Choose the best clustering
        ######################
        self.timer.start("Selection")
        best_clustering_id, all_scores = BestClusteringSelector(clustering_parameters).choose_best(selected_clusterings)
        self.timer.stop("Selection")

        return best_clustering_id, selected_clusterings, not_selected_clusterings, all_scores
